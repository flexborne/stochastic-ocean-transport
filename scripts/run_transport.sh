#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-../configs/transport_demo.nml}
cd "$(dirname "$0")/../fortran"
make FC=${FC:-gfortran}
./transport_solver "$CONFIG"
