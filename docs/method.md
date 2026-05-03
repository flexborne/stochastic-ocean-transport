# Method overview

This repository contains research code for an idealized ensemble simulation of ocean contaminant transport.

## Workflow

1. Generate correlated random perturbations of the wind forcing amplitude.
2. Expand the perturbed forcing in Fourier modes.
3. Compute an ensemble of Stommel stream-function fields.
4. Derive velocity components from the stream function.
5. Run the contaminant transport solver for each velocity realization.
6. Compute ensemble statistics: mean fields, dispersion fields, and spatial correlations.

## Transport model

The contaminant is treated as a passive scalar. Its concentration is governed by a two-dimensional advection-diffusion equation with a short-duration point source. The numerical approach follows a fractional-step method: diffusion is treated by an explicit central-difference scheme, while advection is solved using a TVD Lax-Wendroff-type method with a flux limiter.

## Data assimilation context

The related publication uses the generated ensemble to test data assimilation of concentration measurements. This repository contains the circulation-generation workflow, transport-model entry point, configuration files, and statistical post-processing code.
