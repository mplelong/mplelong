#!/bin/csh
#
# LSF batch script to run an MPI application 
#
#BSUB -P UCSD0019          # project cod
#BSUB -W 12:00             # wall-clock time (hrs:mins)
#BSUB -n 256     	   # number of MPI tasks in job`
#BSUB -R "span[ptile=16]"  # run 16 MPI tasks per nod
#BSUB -J myjob             # job nam
#BSUB -o myjob.3d_test.out # output file name
#BSUB -e myjob.3d_test.err # error file name
#BSUB -q regular           # queue: small, regular


#------------------------------------------------------------------------
#  cd to the run directory on the lustre filesystem for this run
#    $WORK      no backups, not purged, quota 512Go
#    $SCRATCH   no backups, 10TB, purged after 60 days of no use
#------------------------------------------------------------------------
export RUNDIR=/glade/scratch/marinet/resol2/NT
cd $RUNDIR
#run the executable 
mpirun.lsf ./flow.x >& runlog
