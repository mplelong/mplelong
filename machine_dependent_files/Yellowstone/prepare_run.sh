#!/bin/bash
#------------------------------------------------------------------------
# load modules needed for compilation
#------------------------------------------------------------------------
#------------------------------------------------------------------------
module load fftw mkl python all-python-libs


#------------------------------------------------------------------------
# compile and create output directories in my $HOME directory
#------------------------------------------------------------------------
make clean
make flow.x
make outdirs

#------------------------------------------------------------------------
#  create a new directory on the lustre filesystem for this run
#    $WORK      no backups, not purged, quota, 512Go
#    $SCRATCH   no backups,  PURGED after 60 days, 10TB
#------------------------------------------------------------------------
export RUNDIR=/glade/scratch/marinet/resol2/NT
rm -rf $RUNDIR
mkdir $RUNDIR

#------------------------------------------------------------------------
# MOVE flow.x and the RESTART and output directories to the run directory
# COPY the (symbolically linked) input directory to the run directory
#------------------------------------------------------------------------
mv flow.x $RUNDIR/.
mv output $RUNDIR/.
cp -rH input $RUNDIR/.
#cp -rH RESTART $RUNDIR/input/. 
cp -f run.sh Make.inc Makefile $RUNDIR/.

#---------------------------------------------------------------------------------
# cd to $RUNDIR and create the restart files
# this assumes that the 2d .nc file is in input/data_tools and that
# input/data_tools/create_restart_files.py has been set up correctly
# the results will reside in $RUNDIR/flow_solve/input/RESTART and this
# should be consistent with the path specified in input/user_params_module.f90
cd $RUNDIR/input/data_tools
ln -sf ../../../../RESTART/XYZ_resol1.nc .
python interp_V.py
python interp_H.py
python XYZ_to_p1p2_restart.py
