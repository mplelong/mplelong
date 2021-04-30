#!/bin/bash
#
#PBS -N pacific_qsub_test
#PBS -l walltime=00:01:00

#---------------------------------------------------------------------
### Select 4 nodes with 40 CPUs each for a total of 160 MPI processes
#---------------------------------------------------------------------
#PBS -l nodes=4:ppn=12

#---------------------------------------------------------------------
### Send email on abort, begin and end
#---------------------------------------------------------------------
##PBS -m abe
##PBS -M kraig@coast.ucsd.edu

##PBS -j oe
##PBS -o runlog

#-------------------------------------------------------------------------
### Run the executable
#-------------------------------------------------------------------------
module load python/python-3.70-anaconda-no-avx fftw/fftw-3.3.8 intel/parallel_studio_xe_2020.0 nco-4.9.2 lapack/lapack-3.8.0 mpi/openmpi-3.1.0-intel
mpirun date >& KBW/runlog




#----------------------------------------------------------------------------------------
### Run the executable
#
#   mpirun with -np X flow.x argument will run X instances of flow.x on a single compute
#   node using X cores where X <= 40
#
#   how do I run say 160 processes using 4 compute nodes?
#   it seems it should be by specifying  #PBS -l nodes=4:ppn=40
#   and then launching flow.x via mpi  W/O further specifying 160 processes
#----------------------------------------------------------------------------------------

#  by cd'ing to $RUNDIR   
#  flow.x, input and output directories are all visible from compute nodes
#export RUNDIR=/atlantic/kucla/nesting_test_1
#cd $RUNDIR

#which mpirun >& mpirun_version
#hostname >& compute_host

#mpirun -np 40 date  >& runlog
#mpirun -np 40 date  >& output/runlog
#mpirun -np 40 date  >& input/runlog
