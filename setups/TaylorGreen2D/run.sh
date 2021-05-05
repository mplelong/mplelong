#!/bin/bash
#-------------------------------------------------------------------------
# Script to compile flow_solve, create a clean output area, run the codes
# and do some basic data management tasks
#-------------------------------------------------------------------------
echo run.sh now running on "$HOSTNAME"

#------------------------------------------------------------------------
# root directory absolute pathname, able to deal w/ symlinks etc
#------------------------------------------------------------------------
flow_solve_root="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

#------------------------------------------------------------------------
# paths to binaries on local machine
#------------------------------------------------------------------------
PYTHON="/usr/local/bin/python3"
MPIRUN="/usr/local/bin/mpirun"
VIEW="/usr/bin/open -a Preview"

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

#---------------------------------------------------------
# when run is complete, stitch together the output files
# make a simple plot of the cfl history
# 2nd arg for plot_cfl.py is s or hrs or days
#---------------------------------------------------------
 input/data_tools
"$PYTHON" concatenate_results.py "$flow_solve_root"
"$PYTHON" python_scripts/plot_cfl.py "$flow_solve_root"  s    # plot the time scale seconds

cd "$flow_solve_root"
/usr/bin/open -a Preview output/figures/cfl.png
