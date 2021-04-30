#!/bin/bash
#---------------------------------------------------------------------
#   job name 
#---------------------------------------------------------------------
#PBS -N p3_c2_6ipers_data

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
#   alternate syntax:   #PBS -lselect=4:ncpus=40:mpiprocs=40
#---------------------------------------------------------------------
#PBS -lnodes=2:ppn=1

#---------------------------------------------------------------------
### Send email on abort, begin and end
#---------------------------------------------------------------------
#PBS -m abe
#PBS -M kraig@coast.ucsd.edu


#-------------------------------------------------------------------------
###  Specify the directory on /atlantic (120 TB) where the executable is
###  NB both flow.x and the output directory "output" exist in $RUNDIR
###  cd to the top level directory for this run
#-------------------------------------------------------------------------
export RUNDIR=/atlantic/kucla/p3_c2_6ipers
cd $RUNDIR
cp compute_nodes $RUNDIR/input/data_tools
cd $RUNDIR/input/data_tools

#-------------------------------------------------------------------------
### Run the executable sending standard out/error to the file runlog
#-------------------------------------------------------------------------
#  XYZ concats taka about 45 minutes per XYZ block but can run concurrently
#  putting a few processes on each node to ensure large memory is available
#  the final files are 4.8 GB   #PBS -lnodes=4:ppn=2
#mpirun -np 8 -machinefile compute_nodes python concat_XYZ_mpi.py  >& /atlantic/kucla/p3_c2_6ipers/XYZ_log

mpirun -np 2 -machinefile compute_nodes python concat_XZ_mpi.py  >& /atlantic/kucla/p3_c2_6ipers/XZ_log





