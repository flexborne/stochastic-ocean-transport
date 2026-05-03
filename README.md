# Stochastic Ocean Tracer Transport

Numerical framework for ensemble simulations of passive tracer transport in an idealized stochastic ocean circulation.

The project combines a stochastic Stommel-type circulation model with a finite-difference advection-diffusion solver. The workflow generates perturbed stream-function fields, derives velocity components, computes ensemble statistics, and simulates contaminant dispersion from a short-duration point source.

## Scientific background

The model setup follows the idealized numerical experiment described in:

> Kovalets, I., Kim, K. O., Shrubkovsky, O., & Maderich, V. (2022). *Ensemble Data Assimilation of Concentration Measurements Following the Accidental Release of a Contaminant in the Ocean: Method Testing in an Idealized Setting*. Pure and Applied Geophysics. https://doi.org/10.1007/s00024-022-02990-5

The circulation field is generated from a Stommel-type stream-function equation in a rectangular basin. Random correlated perturbations are added to the wind-stress curl to produce an ensemble of velocity-field realizations. The passive contaminant is then transported by the corresponding velocity fields according to a two-dimensional advection-diffusion equation with a short-duration point source.

## Repository structure

```text
configs/      YAML and Fortran namelist configuration files
src/          Python package for circulation generation and statistical analysis
fortran/      Fortran advection-diffusion transport solver
scripts/      Command-line scripts for running workflow stages
data/         Local input/output folders; generated data is ignored by Git
docs/         Method notes and appendix-based mathematical description
```

## Pipeline

```text
1. Generate stochastic circulation ensemble
   Input:  configs/demo.yaml or configs/default.yaml
   Code:   scripts/generate_circulation.py, src/circulation/*
   Output: data/output/grid*.txt
           data/output/deviation*.txt
           data/velocity/<member>_u.txt
           data/velocity/<member>_v.txt

2. Analyze generated ensemble
   Input:  generated grid*.txt files
   Code:   scripts/analyze_velocity.py, src/analysis/*
   Output: mean, dispersion, and correlation text fields under data/output/

3. Run contaminant transport
   Input:  configs/transport_demo.nml or configs/transport.nml
           data/velocity/<member>_u.txt
           data/velocity/<member>_v.txt
   Code:   fortran/transport_solver.f90
   Output: data/output/dist_<member>_<day>.dat
```

## Dependencies

Python scripts use:

```text
numpy
scipy
sympy
matplotlib
PyAstronomy
PyYAML
```

The Fortran solver requires `gfortran` or another Fortran compiler.

## Quick local test

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Generate a small demo circulation ensemble and velocity fields:

```bash
python scripts/generate_circulation.py --config configs/demo.yaml
```

Check that files were created:

```bash
ls data/output
ls data/velocity
```

Expected velocity files include:

```text
data/velocity/0_u.txt
data/velocity/0_v.txt
```

Run statistical analysis:

```bash
python scripts/analyze_velocity.py --config configs/demo.yaml
```

Build and run the Fortran transport solver on the demo configuration:

```bash
cd fortran
make FC=gfortran
./transport_solver ../configs/transport_demo.nml
```

Expected transport output:

```text
data/output/dist_0_10.dat
```

## Full-size run

Generate the full configured ensemble:

```bash
python scripts/generate_circulation.py --config configs/default.yaml
```

Build and run the transport solver:

```bash
cd fortran
make FC=gfortran
./transport_solver ../configs/transport.nml
```

The full configuration is heavier than the demo configuration and may take substantially longer.

## Configuration files

Main YAML configuration files:

```text
configs/demo.yaml       Small grid and small ensemble for local tests
configs/default.yaml    Full experiment configuration
```

Fortran namelist files:

```text
configs/transport_demo.nml   Demo transport run, compatible with configs/demo.yaml
configs/transport.nml        Full transport run, compatible with configs/default.yaml
```

Important configurable parameters include:

- basin size and grid resolution;
- number of Fourier modes;
- stochastic forcing amplitude, standard deviation, and correlation radius;
- ensemble size and random seed;
- contaminant diffusivity;
- release duration and source location;
- simulation period;
- input/output directories.

## Output data

Generated data is written under `data/` and is not committed to the repository:

```text
data/input/      optional initial concentration fields for restart runs
data/velocity/   generated velocity fields used by the Fortran solver
data/output/     stream functions, deviations, statistics, and concentration files
data/reference/  optional reference fields for statistical comparison
```

## Notes

The Fortran solver reads velocity fields generated by the Python circulation workflow. The demo namelist uses member `0`, so run the Python generation step before running `transport_solver`.
