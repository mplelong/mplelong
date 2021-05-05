#!/bin/bash
#-------------------------------------------------------------------------
# Script to compile flow_solve, create a clean output area, run the codes
# and do some basic data management tasks
#-------------------------------------------------------------------------
echo run.sh now running on "$HOSTNAME"

# absolute pathname, able to deal w/ symlinks etc
flow_solve_root="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

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
mpirun -np 4 ./flow.x

#---------------------------------------------------------
# when run is complete, stitch together the output files
# make a simple plot of the cfl history
# 2nd arg for plot_cfl.py is s or hrs or days
#---------------------------------------------------------
cd input/data_tools
python concatenate_results.py "$flow_solve_root/"
python python_scripts/plot_cfl.py "$flow_solve_root/"  s    # plot the time scale seconds

cd "$flow_solve_root"
open -a Preview output/figures/cfl.eps
