#!/bin/bash
#------------------------------------------------------------------------
#  basic work flow for flow_solve on comet.sdsc.edu
# (1) ./prepare_run.sh  (2) cd to the run directory  (3) sbatch run.sh
#------------------------------------------------------------------------



#------------------------------------------------------------------------
# First, load modules needed for compilation
#------------------------------------------------------------------------
module purge

module load intel
module load mpi/openmpi-3.1.0-intel
#  ---- for gcc compilers
#module load mpi

module load fftw
module load lapack

#------------------------------------------------------------------------
#    Specify the run directory on /atlantic 120TB disk space, raid 6
#    redundancy but not backed up, NFS mounted on pacific via infiniband
#    physically located on a JBOD everything connected via infiniband
#------------------------------------------------------------------------
export RUNDIR=/atlantic/kucla/p3_c2_6ipers
rm -rf $RUNDIR
mkdir $RUNDIR  


#------------------------------------------------------------------------
# compile and create output directories in my $HOME/... directory
#------------------------------------------------------------------------
cd $HOME/KBW/flow_solve
make clean
make flow.x
make outdirs


#------------------------------------------------------------------------
# MOVE flow.x and output directories to the run directory
# COPY the (symbolically linked) input directory to the run directory
#------------------------------------------------------------------------
mv flow.x $RUNDIR/.
mv output $RUNDIR/.
cp $HOME/KBW/nested_run/4x10/wind_speed_direction input/.
cp -rH input $RUNDIR/.
cp -f run.sh prepare_run.sh compute_nodes process_data.sh test.sh Make.inc Makefile $RUNDIR/.


#---------------------------------------------------------------------------------
# run a script to create or mv restart files to the appropriate location
# if necessary, or set link to saved output from previous run
#----------------------------------------------------------------------------------
ln -s $HOME/KBW/nested_run/8x20/RESTART $RUNDIR/RESTART
ln -s $HOME/KBW/nested_run/8x20/BC $RUNDIR/BC

