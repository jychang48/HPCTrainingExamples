#!/bin/bash

#SBATCH -J vector_addition_hipify
#SBATCH -N 1
#SBATCH -t 5

srun -N1 -n1 ./vector_addition
