#!/bin/bash

module purge
module load netcdf-c

# use the compiler used to build netcdf-c
CC=`nc-config --cc`

rm -rf netcdf-c_test
mkdir netcdf-c_test
cd netcdf-c_test
wget https://www.unidata.ucar.edu/software/netcdf/examples/programs/pres_temp_4D_wr.c
wget https://www.unidata.ucar.edu/software/netcdf/examples/programs/pres_temp_4D_rd.c
$CC pres_temp_4D_wr.c -lnetcdf -L${NETCDF_C_ROOT}/lib -o pres_temp_4D_wr
$CC pres_temp_4D_rd.c -lnetcdf -L${NETCDF_C_ROOT}/lib -o pres_temp_4D_rd
./pres_temp_4D_wr
./pres_temp_4D_rd
cd ..
rm -rf netcdf-c_test



