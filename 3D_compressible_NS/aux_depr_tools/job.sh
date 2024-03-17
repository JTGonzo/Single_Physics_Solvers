#!/bin/bash
#PBS -l walltime=24:00:00,select=1:ncpus=1:mpiprocs=1:mem=16gb
#PBS -N matBub3D_m1
#PBS -A ex-rjaiman-1
#PBS -m abe
#PBS -M suraj.kashyap@ubc.ca
#PBS -o output.txt
#PBS -e error_output.txt
#
# ##################
# PBS example script
# ##################
module load gcc/5.4.0 matlab/R2019b

cd /scratch/ex-rjaiman-1/suraj_cases/Cavitation/matlab3D/m1

matlab -nodisplay -r "main"
~                                                                                                                                                                                                           
~                                                                                                                                                                                                      
~                                                                                           
