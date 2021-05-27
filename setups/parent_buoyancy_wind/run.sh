#!/bin/bash
#-------------------------------------------------------------------------
# Script to compile flow_solve, create a clean output area, run the codes
# and do some basic data management tasks
#-------------------------------------------------------------------------
echo run.sh now running on "$HOSTNAME"

source set_env_variables.sh

do_run=True
do_movie_frames=False
make_restart_files=False

#------------------------------------------------------------
# prepare things for the parent_buoyancy_wind run
#------------------------------------------------------------
if( "$do_run" ) then
	rm -f input
	ln -s setups/parent_buoyancy_wind input
	make clean
	make flow.x
	make outdirs

	#---------------------------------------------------
	# do the run (it's set up for a 2x2 processor grid)
	#---------------------------------------------------
	"$MPIRUN" -np 4 ./flow.x
fi

#---------------------------------------------------------------
# when run is complete, stitch together the output files
# make a simple plot of the cfl history
# 2nd arg for plot_cfl.py is s, hrs or days for plot time scale
#---------------------------------------------------------------
cd input/data_tools
	"$PYTHON" concatenate_results.py "$FLOW_SOLVE_ROOT"
	"$PYTHON" python_scripts/plot_cfl.py "$FLOW_SOLVE_ROOT"  hrs  # plot time scale hrs

cd "$FLOW_SOLVE_ROOT"
	open -a Preview output/figures/cfl.pdf

	ncrcat -O output/slices/2D/XY* output/slices/2D/XY.nc
	rm -f output/slices/2D/XY_*.nc

	ncrcat -O output/slices/2D/YZ* output/slices/2D/YZ.nc
	rm -f output/slices/2D/YZ_*.nc

	ncrcat -O output/slices/2D/XZ* output/slices/2D/XZ.nc
	rm -f output/slices/2D/XZ_*.nc
	"$NCVIEW" output/slices/2D/YZ.nc &

if( "$do_movie_frames" ) then
	cd input/data_tools/python_scripts
	"$PYTHON" make_movie_frames.py "$FLOW_SOLVE_ROOT"

	cd "$FLOW_SOLVE_ROOT"
	ffmpeg -framerate 10 -pattern_type glob -i 'output/movie_frames/U*.png' -c:v libx264 -vf "fps=10,format=yuv422p" output/XZ_U_nested.avi
fi

#---------------------------------------------------------------
# make a 2x2 set of restart files from global file
# output/slices/3D/XYZ_000768.nc
#---------------------------------------------------------------

if( "$make_restart_files" ) then
 	"$PYTHON" "$FLOW_SOLVE_ROOT"/input/data_tools/python_scripts/XYZ_to_p1p2_restart.py "$FLOW_SOLVE_ROOT" 2 2 XYZ_000768.nc
fi



