#!/bin/bash

module load gcc/13

REPO_DIR="$(dirname "$(dirname "$(readlink -fm "$0")")")"
cd ${REPO_DIR}/Pragma_Examples/OpenMP/Fortran/vecadd

make
./vecadd
make clean
