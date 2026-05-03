# Analytical solution for the stream function

The circulation ensemble is generated using an idealized linear Stommel model for wind-driven circulation in a rectangular basin. The stream function satisfies a stationary linear equation with no-flux boundary conditions on the lateral boundaries.

For a constant wind forcing amplitude, the solution is written as the sum of a homogeneous and a particular part. For a perturbed forcing amplitude, the right-hand side is expanded into a Fourier series. Each Fourier mode contributes a separable component to the stream function field.

In the numerical implementation, the infinite expansion is truncated to a finite number of modes. The paper appendix reports that 15 modes were sufficient because higher-order terms had negligible influence on the stream-function field.

The velocity components are recovered from the stream function as

```text
u = dpsi/dy
v = -dpsi/dx
```

These velocity fields can then be used by the advection-diffusion solver for passive contaminant transport.
