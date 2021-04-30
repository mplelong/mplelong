#!/bin/bash
#------------------------------------------------------------------------
#  basic work flow for flow_solve on comet.sdsc.edu
# (1) ./prepare_run.sh  (2) cd to the run directory  (3) sbatch run.sh
#------------------------------------------------------------------------



#------------------------------------------------------------------------
# First, load modules needed for compilation
#------------------------------------------------------------------------
#module purge
#module load gnutools
#module load intel mvapich2_ib
##module load pgi mvapich2_ib
#module load netcdf hdf5 fftw atlas

#------------------------------------------------------------------------
#    Specify the run directory on the lustre filesystem for this run
#    $PROJECTS      no backups, not purged, quota
#    $SCRATCH       no backups,  maybe PURGED after 10 days,  much larger
#------------------------------------------------------------------------
export RUNDIR=$PROJECTS/LockExchange/9.0_run_2
rm -rf $RUNDIR
mkdir $RUNDIR  


#------------------------------------------------------------------------
# compile and create output directories in my $HOME/... directory
#------------------------------------------------------------------------
cd $HOME/flow_solve/flow_solve_GIT_REPOSITORY
make clean
make flow.x
make outdirs


#------------------------------------------------------------------------
# MOVE flow.x and output directories to the run directory
# COPY the (symbolically linked) input directory to the run directory
#------------------------------------------------------------------------
mv flow.x $RUNDIR/.
mv output $RUNDIR/.
cp -rH input $RUNDIR/.
cp -f run.sh prepare_run.sh process_data.sh Make.inc Makefile $RUNDIR/.


#---------------------------------------------------------------------------------
# run a script to create or mv restart files to the appropriate location
# if necessary, or set link to saved output from previous run
#----------------------------------------------------------------------------------
ln -s $PROJECTS/LockExchange/9.0_run_1/output/3D/ $RUNDIR/RESTART

