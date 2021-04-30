#!/bin/bash
#
#SBATCH -J basic_processing         # job name
#SBATCH -o msgs%j                   # output and error file name (%j expands to jobID)
#SBATCH --nodes=1                   # total number of NODES  (24 cores/node)
#SBATCH --ntasks-per-node=1         # total number of mpi tasks requested
#SBATCH --partition=large-shared          # queue (partition) -- compute,shared,large-shared
#SBATCH --mem=64G
#SBATCH --export=ALL
#SBATCH -t 02:00:00                 # run time (hh:mm:ss)  168hrs=1week
##SBATCH --mail-user=kraig@coast.ucsd.edu
##SBATCH --mail-type=begin           # email me when the job starts
##SBATCH --mail-type=end             # email me when the job finishes

##  compute nodes have 128 GB memory  (newton thin have only 64 GB)
##  there are 4 large memory nodes with 1455 GB memory
## For example, on the "large-shared" partition, the following job requesting 16 cores 
## and 445 GB of memory (about 31.3% of 1455 GB of one node's available memory) 
## for 1 hour will be charged 20 SUs:
## 455/1455(memory) * 64(cores) * 1(duration) ~= 20

##SBATCH --ntasks=16
##SBATCH --mem=455G
##SBATCH --partition=large-shared

#--------------------------------------------------------------------------
# commands for basic data processing
# NB no specification of the number of  processes here
#   (1) configure the perl and python scripts in input/data_tools
#   (2) from the flow_solve directory:  sbatch process_data.sh
#--------------------------------------------------------------------------

rm msgs*
cd input/data_tools
chmod a+x *.pl

#rm -f cfl_log; python plot_cfl.py < /dev/null >& cfl_log     # do this routinely

#---------------------------------------------------------------
# concatenation of distributed data files only
#---------------------------------------------------------------
#ibrun python concat_YZ_mpi.py >& YZ_concat_log             # run the python script using ntasks-per-node  
#ibrun python concat_XYZ_mpi.py >& XYZ_concat_log           # run the python script using ntasks-per-node
#ibrun python concat_XY_low_mpi.py >& XY_low_concat_log     # run the python script using ntasks-per-node   
#ibrun python concat_XY_high_mpi.py >& XY_high_concat_log   # run the python script using ntasks-per-node    
#  XY: don't let the processors compete w/ one another for access to the same files  NECESSARY???

#---------------------------------------------------------------
# optional concatenation and plotting
#---------------------------------------------------------------
#rm -f YZ_log; ibrun python plot_YZ_slice.py < /dev/null >& YZ_log
#rm -f XY_log; ibrun python plot_XY_slice.py < /dev/null >& XY_log
#rm -f XY_high_log; ibrun python plot_XY_high_slice.py < /dev/null >& XY_high_log

#  ~ 6 minutes to do the 3d joins,  with 64GB    ~8 minutes total    5:15 w/ create 3D
rm -f sorted_log energy.dat; ibrun python sorted_profile.py < /dev/null >& sorted_log



