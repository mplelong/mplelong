#!/bin/bash
#-------------------------------------------------------------------------
# Script to compile flow_solve, create a clean output area, run the codes
# and do some basic data management tasks
#-------------------------------------------------------------------------
echo run.sh now running on lacosta

#------------------
# prepare things
#------------------
make clean
make flow.x
make outdirs

#------------------
# do the run
#------------------
mpirun -np 4 ./flow.x

#---------------------------------------------------------
# when run is complete, stitch together the output files
# make a simple plot of the cfl history
# 2nd arg for plot_cfl.py is s or hrs or days
#---------------------------------------------------------
cd input/data_tools
python concatenate_results.py
python python_scripts/plot_cfl.py /Users/kraig/flow_solve_BC/ s    # plot the time scale seconds

cd /Users/kraig/flow_solve_BC/
open -a Preview output/figures/cfl.eps
