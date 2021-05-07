#!/bin/bash
#-------------------------------------------------------------------------
# Script to compile flow_solve, create a clean output area, run the codes
# and do some basic data management tasks
#-------------------------------------------------------------------------
echo run.sh now running on "$HOSTNAME"

source set_env_variables.sh

#-------------------------------------------------
# prepare things for Taylor-Green vortex test run
#-------------------------------------------------
rm -f input
ln -s setups/TaylorGreen2D input
make clean
make flow.x
make outdirs

#---------------------------------------------------
# do the run (it's set up for a 2x2 processor grid)
#---------------------------------------------------
"$MPIRUN" -np 4 ./flow.x

#---------------------------------------------------------------
# when run is complete, stitch together the output files
# make a simple plot of the cfl history
# 2nd arg for plot_cfl.py is s, hrs or days for plot time scale
#---------------------------------------------------------------
cd input/data_tools
"$PYTHON" concatenate_results.py "$FLOW_SOLVE_ROOT"
"$PYTHON" python_scripts/plot_cfl.py "$FLOW_SOLVE_ROOT"  s    # plot the time scale seconds

cd "$FLOW_SOLVE_ROOT"
open -a Preview output/figures/cfl.pdf
