#!/bin/bash

module purge
module load netcdf-c
module load netcdf-fortran
module load openmpi

# use the compiler used to build netcdf-fortran
FC=`nf-config --fc`

rm -rf netcdf-fortran_test
mkdir netcdf-fortran_test
cd netcdf-fortran_test
git clone https://github.com/Unidata/netcdf-fortran.git
$FC  ./netcdf-fortran/examples/F90/simple_xy_par_wr.F90 -o simple_xy_par_wf  -I${NETCDF_FC_ROOT}/include -L${NETCDF_FC_ROOT}/lib -lnetcdff -L${PNETCDF_ROOT}/lib -lpnetcdf
mpirun -n 4 ./simple_xy_par_wf
ncdump simple_xy_par.nc
cd ..
rm -rf netcdf-fortran_test



