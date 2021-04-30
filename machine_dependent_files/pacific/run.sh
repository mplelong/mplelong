#!/bin/bash
#---------------------------------------------------------------------
#   job name 
#---------------------------------------------------------------------
#PBS -N pacific_nesting_test_4

#---------------------------------------------------------------------
#   run time (hh:mm:ss) 
#---------------------------------------------------------------------
#PBS -l walltime=48:00:00

#---------------------------------------------------------------------
#   name of queue 
#---------------------------------------------------------------------				
#PBS -q workq

#---------------------------------------------------------------------
#   Select 4 nodes with 40 CPUs each for a total of 160 MPI processes
#    alternate syntax:   #PBS -lselect=4:ncpus=40:mpiprocs=40
#---------------------------------------------------------------------
#PBS -lnodes=4:ppn=40

#---------------------------------------------------------------------
### Send email on abort, begin and end
#---------------------------------------------------------------------
#PBS -m abe
#PBS -M kraig@coast.ucsd.edu

#-------------------------------------------------------------------------
# load modules 
#   NB netcdf-fortran installed by me and so explicit links to
#      libraries for compilation specified in Make.inc. However,
#      since netcdf is not under module control, the shared library
#      needed is not in my path and so needs to be added explicitly
#-------------------------------------------------------------------------
module purge
#  --- for intel mpif90
#module load intel
#module load mpi/openmpi-3.1.0-intel
#module load fftw
#module load lapack

module load intel/parallel_studio_xe_2020.0
module load hdf5/hdf5-1.10.6_intel
module load netcdf-c-4.7.3_intel_mpicc_mpif90
module load netcdf-fortran-4.5.2_intel_mpicc_mpif90
module load mpi/openmpi-3.1.0-intel
module load fftw
module load lapack


#-------------------------------------------------------------------------
###  Specify the directory on /atlantic (120 TB) where the executable is
###  NB both flow.x and the output directory "output" exist in $RUNDIR
###  cd to the top level directory for this run
#-------------------------------------------------------------------------
export RUNDIR=/atlantic/kucla/p3_c2_6ipers
cd $RUNDIR

#-------------------------------------------------------------------------
### Run the executable sending standard out/error to the file runlog
#-------------------------------------------------------------------------
mpirun -np 160 -machinefile compute_nodes ./flow.x >& runlog
date >> runlog




