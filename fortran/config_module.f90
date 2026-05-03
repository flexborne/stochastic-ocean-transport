module config_module
  implicit none

  integer :: nx = 64
  integer :: ny = 64
  integer :: time_steps_per_day = 16
  integer :: member_id = 0
  integer :: day_start = 0
  integer :: day_end = 10
  integer :: source_i = 10
  integer :: source_j = 35

  double precision :: length_x = 1.0d7
  double precision :: length_y = 6.249d6
  double precision :: diffusivity = 1.0d3
  double precision :: seconds_per_day = 86400.0d0
  double precision :: release_duration_days = 10.0d0
  double precision :: release_rate = 1.0d0

  character(len=256) :: velocity_dir = '../data/velocity'
  character(len=256) :: input_dir = '../data/input'
  character(len=256) :: output_dir = '../data/output'

contains

  subroutine read_transport_config(path)
    character(len=*), intent(in) :: path
    integer :: unit_id, ios

    namelist /transport_config/ nx, ny, length_x, length_y, diffusivity, &
      seconds_per_day, time_steps_per_day, member_id, day_start, day_end, &
      source_i, source_j, release_duration_days, release_rate, &
      velocity_dir, input_dir, output_dir

    open(newunit=unit_id, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Warning: could not open config file: ', trim(path)
      print *, 'Using built-in defaults.'
      return
    end if
    read(unit_id, nml=transport_config)
    close(unit_id)
  end subroutine read_transport_config

end module config_module
