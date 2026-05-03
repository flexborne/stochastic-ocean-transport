program transport_solver
  use config_module
  implicit none

  character(len=256) :: config_path

  if (command_argument_count() >= 1) then
    call get_command_argument(1, config_path)
  else
    config_path = 'transport.nml'
  end if

  call read_transport_config(trim(config_path))
  call print_configuration()
  call run_transport()

contains

  subroutine print_configuration()
    print *, 'Transport solver configuration:'
    print *, '  nx, ny = ', nx, ny
    print *, '  member_id = ', member_id
    print *, '  day_start, day_end = ', day_start, day_end
    print *, '  length_x, length_y = ', length_x, length_y
    print *, '  diffusivity = ', diffusivity
    print *, '  source_i, source_j = ', source_i, source_j
    print *, '  velocity_dir = ', trim(velocity_dir)
    print *, '  input_dir = ', trim(input_dir)
    print *, '  output_dir = ', trim(output_dir)
  end subroutine print_configuration

  subroutine run_transport()
    integer :: i, j, k, n_steps, current_layer, diff_layer, next_layer
    integer :: write_day_index, unit_id, ios
    integer, dimension(9) :: checkpoints
    double precision :: dx, dy, dt, elapsed_day
    double precision, allocatable :: c(:, :, :), c_star(:, :), c_initial(:, :)
    double precision, allocatable :: u(:, :), v(:, :), array_out(:, :, :)
    character(len=512) :: u_path, v_path, output_path, input_path

    checkpoints = (/100, 200, 300, 400, 500, 600, 700, 800, 900/)

    if (day_end <= day_start) then
      print *, 'Error: day_end must be greater than day_start.'
      stop 1
    end if

    dx = length_x / dble(nx - 1)
    dy = length_y / dble(ny - 1)
    n_steps = (day_end - day_start) * time_steps_per_day + 1
    dt = seconds_per_day * dble(day_end - day_start) / dble(n_steps - 1)

    allocate(c(nx, ny, 3))
    allocate(c_star(nx, ny))
    allocate(c_initial(nx, ny))
    allocate(u(nx, ny))
    allocate(v(nx, ny))
    allocate(array_out(1, nx, ny))

    c = 0.0d0
    c_star = 0.0d0
    c_initial = 0.0d0
    u = 0.0d0
    v = 0.0d0

    write(u_path, '(A, "/", I0, "_u.txt")') trim(velocity_dir), member_id
    write(v_path, '(A, "/", I0, "_v.txt")') trim(velocity_dir), member_id
    call read_velocity_file(trim(u_path), u)
    call read_velocity_file(trim(v_path), v)
    call apply_boundary_velocity(u, v)

    if (day_start > 0) then
      write(input_path, '(A, "/dist_", I0, "_", I0, ".dat")') trim(input_dir), member_id, day_start
      call read_concentration(trim(input_path), c(:, :, 1))
      c_initial = c(:, :, 1)
    else
      c(:, :, :) = 0.0d0
      c_initial = 0.0d0
    end if

    print *, 'Running transport for ', n_steps - 1, ' time steps; dt = ', dt

    do k = 1, n_steps - 1
      current_layer = mod(k + 2, 3) + 1
      diff_layer = current_layer
      next_layer = mod(k, 3) + 1
      elapsed_day = dble(day_start) + dt * dble(k) / seconds_per_day

      if (day_start == 0 .and. dt * dble(k) <= release_duration_days * seconds_per_day) then
        if (source_i >= 1 .and. source_i <= nx .and. source_j >= 1 .and. source_j <= ny) then
          c(source_i, source_j, current_layer) = c(source_i, source_j, current_layer) + release_rate * dt
        end if
      end if

      do i = 2, nx - 1
        do j = 2, ny - 1
          c_star(i, j) = diffusion_step(i, j, current_layer, c, dx, dy, dt)
        end do
      end do

      c(:, :, diff_layer) = 0.0d0
      do i = 2, nx - 1
        do j = 2, ny - 1
          c(i, j, diff_layer) = c_star(i, j)
        end do
      end do

      c(:, :, next_layer) = 0.0d0
      do i = 2, nx - 1
        do j = 2, ny - 1
          c(i, j, next_layer) = advection_step(i, j, diff_layer, c, u, v, dx, dy, dt)
        end do
      end do

      do write_day_index = 1, size(checkpoints)
        if (abs(elapsed_day - dble(checkpoints(write_day_index))) < 0.5d0 / dble(time_steps_per_day)) then
          if (checkpoints(write_day_index) > day_start .and. checkpoints(write_day_index) < day_end) then
            write(output_path, '(A, "/dist_", I0, "_", I0, ".dat")') trim(output_dir), member_id, checkpoints(write_day_index)
            call write_concentration(trim(output_path), c(:, :, next_layer), array_out)
          end if
        end if
      end do

      if (mod(k, max(1, (n_steps - 1) / 10)) == 0) then
        print *, '  completed step ', k, ' / ', n_steps - 1
      end if
    end do

    write(output_path, '(A, "/dist_", I0, "_", I0, ".dat")') trim(output_dir), member_id, day_end
    call write_concentration(trim(output_path), c(:, :, mod(n_steps - 1, 3) + 1), array_out)
    call print_norms(c(:, :, mod(n_steps - 1, 3) + 1), c_initial, dx, dy)
    print *, 'Transport output written to ', trim(output_path)

    deallocate(c, c_star, c_initial, u, v, array_out)
  end subroutine run_transport

  double precision function limiter(s)
    double precision, intent(in) :: s
    limiter = (abs(s) + s) / (abs(s) + 1.0d0)
    if (s < 0.0d0) limiter = 0.0d0
  end function limiter

  double precision function flux(i, j, layer, direction, c, u, v, dx, dy, dt)
    integer, intent(in) :: i, j, layer, direction
    double precision, intent(in) :: c(:, :, :), u(:, :), v(:, :), dx, dy, dt
    double precision :: dp, s, velocity_value, spacing
    integer :: shift, i0, i1, j0, j1

    if (direction == 1) then
      dp = c(i + 1, j, layer) - c(i, j, layer)
      velocity_value = u(i, j)
      spacing = dx
      if (abs(dp) <= 1.0d-50) then
        s = 0.0d0
      else
        shift = int(sign(1.0d0, -velocity_value))
        i1 = max(1, min(nx, i + 1 + shift))
        i0 = max(1, min(nx, i + shift))
        s = (c(i1, j, layer) - c(i0, j, layer)) / dp
      end if
      flux = velocity_value * (c(i + 1, j, layer) + c(i, j, layer)) / 2.0d0
      flux = flux - dp * ((1.0d0 - limiter(s)) * abs(velocity_value) + velocity_value * velocity_value * dt * limiter(s) / spacing) / 2.0d0
    else
      dp = c(i, j + 1, layer) - c(i, j, layer)
      velocity_value = v(i, j)
      spacing = dy
      if (abs(dp) <= 1.0d-50) then
        s = 0.0d0
      else
        shift = int(sign(1.0d0, -velocity_value))
        j1 = max(1, min(ny, j + 1 + shift))
        j0 = max(1, min(ny, j + shift))
        s = (c(i, j1, layer) - c(i, j0, layer)) / dp
      end if
      flux = velocity_value * (c(i, j + 1, layer) + c(i, j, layer)) / 2.0d0
      flux = flux - dp * ((1.0d0 - limiter(s)) * abs(velocity_value) + velocity_value * velocity_value * dt * limiter(s) / spacing) / 2.0d0
    end if

    if (abs(flux) < 1.0d-70) flux = 0.0d0
  end function flux

  double precision function advection_step(i, j, layer, c, u, v, dx, dy, dt)
    integer, intent(in) :: i, j, layer
    double precision, intent(in) :: c(:, :, :), u(:, :), v(:, :), dx, dy, dt
    advection_step = c(i, j, layer) &
      - dt * (flux(i, j, layer, 1, c, u, v, dx, dy, dt) - flux(i - 1, j, layer, 1, c, u, v, dx, dy, dt)) / dx &
      - dt * (flux(i, j, layer, 2, c, u, v, dx, dy, dt) - flux(i, j - 1, layer, 2, c, u, v, dx, dy, dt)) / dy
  end function advection_step

  double precision function diffusion_step(i, j, layer, c, dx, dy, dt)
    integer, intent(in) :: i, j, layer
    double precision, intent(in) :: c(:, :, :), dx, dy, dt
    diffusion_step = c(i, j, layer)
    if (i > 1 .and. i < nx) then
      diffusion_step = diffusion_step + diffusivity * (c(i + 1, j, layer) - 2.0d0 * c(i, j, layer) + c(i - 1, j, layer)) * dt / (dx * dx)
    end if
    if (j > 1 .and. j < ny) then
      diffusion_step = diffusion_step + diffusivity * (c(i, j + 1, layer) - 2.0d0 * c(i, j, layer) + c(i, j - 1, layer)) * dt / (dy * dy)
    end if
    if (j == 2 .or. j == ny - 1) then
      diffusion_step = diffusion_step + diffusivity * c(i, j, layer) * dt / (dy * dy)
    end if
    if (i == 2 .or. i == nx - 1) then
      diffusion_step = diffusion_step + diffusivity * c(i, j, layer) * dt / (dx * dx)
    end if
  end function diffusion_step

  subroutine apply_boundary_velocity(u, v)
    double precision, intent(inout) :: u(:, :), v(:, :)
    integer :: i
    do i = 3, nx - 2
      if (v(i, 2) < 0.0d0) v(i, 2) = 0.0d0
      if (v(i, ny - 1) > 0.0d0) v(i, ny - 1) = 0.0d0
      if (u(2, i) < 0.0d0) u(2, i) = 0.0d0
      if (u(nx - 1, i) > 0.0d0) u(nx - 1, i) = 0.0d0
    end do
    u(2, 2) = 0.0d0
    u(nx - 1, ny - 1) = 0.0d0
    v(2, ny - 1) = 0.0d0
    v(nx - 1, 2) = 0.0d0
  end subroutine apply_boundary_velocity

  subroutine read_velocity_file(path, field)
    character(len=*), intent(in) :: path
    double precision, intent(out) :: field(:, :)
    integer :: unit_id, ios, row
    open(newunit=unit_id, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: could not open velocity file ', trim(path)
      stop 2
    end if
    do row = 1, ny
      read(unit_id, *, iostat=ios) field(:, row)
      if (ios /= 0) then
        print *, 'Error: failed to read row ', row, ' from ', trim(path)
        stop 3
      end if
    end do
    close(unit_id)
  end subroutine read_velocity_file

  subroutine read_concentration(path, concentration)
    character(len=*), intent(in) :: path
    double precision, intent(out) :: concentration(:, :)
    real, allocatable :: tmp(:, :, :)
    integer :: unit_id, km, im, jm, ios
    open(newunit=unit_id, file=trim(path), status='old', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: could not open restart concentration file ', trim(path)
      stop 4
    end if
    read(unit_id) km, im, jm
    if (im /= nx .or. jm /= ny) then
      print *, 'Error: restart grid does not match config.'
      stop 5
    end if
    allocate(tmp(km, im, jm))
    read(unit_id) tmp
    close(unit_id)
    concentration(:, :) = dble(tmp(1, :, :))
    deallocate(tmp)
  end subroutine read_concentration

  subroutine write_concentration(path, concentration, array_out)
    character(len=*), intent(in) :: path
    double precision, intent(in) :: concentration(:, :)
    double precision, intent(inout) :: array_out(:, :, :)
    integer :: unit_id, ios
    call execute_command_line('mkdir -p ' // trim(output_dir), wait=.true.)
    array_out(1, :, :) = concentration(:, :)
    open(newunit=unit_id, file=trim(path), form='unformatted', action='write', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: could not open output file ', trim(path)
      stop 6
    end if
    write(unit_id) 1, nx, ny
    write(unit_id) array_out(1, :, :)
    close(unit_id)
  end subroutine write_concentration

  subroutine print_norms(concentration, initial, dx, dy)
    double precision, intent(in) :: concentration(:, :), initial(:, :), dx, dy
    double precision :: min_c, max_c, max_abs_delta, l1_delta, l2_delta, mass
    integer :: i, j
    min_c = huge(1.0d0)
    max_c = -huge(1.0d0)
    max_abs_delta = 0.0d0
    l1_delta = 0.0d0
    l2_delta = 0.0d0
    mass = 0.0d0
    do i = 2, nx - 1
      do j = 2, ny - 1
        min_c = min(min_c, concentration(i, j))
        max_c = max(max_c, concentration(i, j))
        max_abs_delta = max(max_abs_delta, abs(concentration(i, j) - initial(i, j)))
        l1_delta = l1_delta + dx * dy * abs(concentration(i, j) - initial(i, j))
        l2_delta = l2_delta + ((concentration(i, j) - initial(i, j)) * dx) ** 2
      end do
    end do
    do i = 1, nx
      do j = 1, ny
        mass = mass + concentration(i, j)
      end do
    end do
    l2_delta = sqrt(l2_delta)
    print *, 'Diagnostics:'
    print *, '  min, max = ', min_c, max_c
    print *, '  max_abs_delta = ', max_abs_delta
    print *, '  l1_delta = ', l1_delta
    print *, '  l2_delta = ', l2_delta
    print *, '  mass = ', mass
  end subroutine print_norms

end program transport_solver
