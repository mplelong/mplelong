#!/bin/bash
#
#SBATCH -J 9.0_run_2                # job name
#SBATCH -o msgs%j                   # output and error file name (%j expands to jobID)
#SBATCH --nodes=48                  # total number of NODES  (24 cores/node)
#SBATCH --ntasks-per-node=24        # total number of mpi tasks requested
#SBATCH --partition=compute         # queue (partition) -- compute,
#SBATCH --export=ALL
#SBATCH -t 07:30:00                 # run time (hh:mm:ss)  168hrs=1week
#SBATCH --mail-user=kraig@coast.ucsd.edu
#SBATCH --mail-type=begin           # email me when the job starts
#SBATCH --mail-type=end             # email me when the job finishes
######SBATCH --qos=oneweek
#-----------------------------------------------------------------------------------------------
# Kraig,
# I have enabled you for long runs on Comet. You can add the following line to your script:
#        #SBATCH --qos=oneweek
# and then just request the longer run time. 
# Let me know if there are any issues with requesting this longer time.
#  -Mahidhar
#-----------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------
# run the simulation, NB no specification of the number of  processes here
#--------------------------------------------------------------------------
ibrun ./flow.x >& runlog    # run the MPI executable named flow.x
date >> runlog
