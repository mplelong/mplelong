#!
"""
Routine to organize the concatenation of the distributed netcdf data files. 
Concatenated files are written to the output/slices/2D and output/slices/2D directories.

In a multi processor flow_solve run, the data making up a single 2d plane or 3d block is 
saved in output/2D and output/3D as separate netcdf files for each processor. This script
glues these files together, creating separate spatially global files for each time slice.
Once this has been done, the output data can be regarded as residing in output/slices
and the original files in output/2D and output3D can be considered redundant and deleted.

Obviously, check things out carefully before deleting.

This script is executed on a single processor but issues execution commands to
worker routines (e.g. concat_XZ_mpi.py) via mpirun that execute the work using 
multiple processors. The processors divide up the time slices between themselves
and sync up befor returning control to this script.
"""
#---------------------------------------------------------------------------------------
#  import and name the various modules I'll use
#---------------------------------------------------------------------------------------
import os,math,sys         
import numpy as np
from python_scripts.data_processing_utilities import parse_problem_params

#--------------------------------------
# script control params
#--------------------------------------
do_XY=True  ; delete_XY=False
do_XZ=True  ; delete_XZ=False
do_YZ=True  ; delete_YZ=False
do_XYZ=True ; delete_XYZ=False

#------------------------------------------
# make a directory to store file locations
#------------------------------------------
if not os.path.exists('paths'):
	cmd = 'mkdir paths'
	os.system(cmd)



#--------------------------------------------------------------------------- 
# flow_solve root directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
root_dir = sys.argv[1]
#  add / for convenience
root_dir = root_dir + '/'

#------------------------------------------------------------------------
# path to python w/ mpi4py and corresponding mpirun
#------------------------------------------------------------------------
PYTHON="$PYTHON"    # get env variable
MPIRUN="$MPIRUN"    # get env variable



#--------------------------------------------------------------------------- 
# get the processor decomposition used by flow_solve p1,p2 
#---------------------------------------------------------------------------
params = parse_problem_params(root_dir)  

[runlabel,restart_flag,do_second_scalar,AB_order,p1,p2,nx,ny,nz,dt,t0,tf,Lx,Ly,Lz, \
x_periodic,y_periodic,z_periodic,s1_def,s2_def,user_forcing_flag,rho0,g,f0,nu,     \
kappa1,kappa2,high_order_flag,p,T_diff] = params



print("---------------------------------------------------------------------------------")
print (" Running concatenate_results.py XY, XZ, YZ, XYZ: ",do_XY,do_XZ,do_YZ,do_XYZ) 
print (" ...  delete distributed files  XY, XZ, YZ, XYZ: ",delete_XY,delete_XZ,delete_YZ,delete_XYZ)
print("---------------------------------------------------------------------------------")



#------------------------------------------------------------------------------ 
# Task 1: concatenate the XZ plane nc files together
# creating a separate file for each time slice.
# Results stored in output/slices/2D.
# e.g.
# Saved filenames: fnrs_iproc-jproc.nc
# XZ_000-000.nc	XZ_000-001.nc	XZ_001-000.nc	XZ_001-001.nc
#-------------------------------------------------------------------------------
if( do_XZ ):
	fnrs = 'XZ'                             # filename root string
	start_slice = 0                         # starting time slice  i.e. slices = np.arange(start_slice,end_slice,inc)
	end_slice = 51                          # ending time slice    
	inc = 1                                 # time slice increment
	numprocs = 4                            # number of processors to be used by concat_XZ_mpi.py	
	
	cmd_base = MPIRUN + ' -np ' + str(numprocs) + ' ' + PYTHON + ' python_scripts/concat_XZ_mpi.py '
	command = cmd_base + root_dir + ' ' + str(p1) + ' ' + str(p2) + ' ' + str(start_slice) + ' ' + str(end_slice) + ' ' + str(inc) + ' ' + fnrs
	os.system(command)
	if( delete_XZ ):
		command = 'rm -f ' + root_dir + 'output/2D/' + fnrs + '_*.nc'
		os.system(command)
	command = 'ls -lh ' + root_dir + 'output/slices/2D/' + fnrs + '* > paths/concatenated_' + fnrs + '_files'
	os.system(command)


#---------------------------------------------------- 
# Task 2: concatenate the YZ plane nc files together
# creating a separate file for each time slice.
# Results stored in output/slices/2D.
# e.g.
# Saved filenames: fnrs_iproc-jproc.nc
# YZ_001-000.nc	YZ_001-001.nc
#----------------------------------------------------
if( do_YZ ):
	fnrs = 'YZ'                             # filename root string
	start_slice = 0                         # starting time slice  i.e. slices = np.arange(start_slice,end_slice,inc)
	end_slice = 51                          # ending time slice    
	inc = 1                                 # time slice increment
	iproc = 1                               # iproc for the YZ plane (depends on x value of YZ plane saved)
	numprocs = 4                            # number of processors to be used by concat_YZ_mpi.py	
	
	cmd_base = MPIRUN + ' -np ' + str(numprocs) + ' ' + PYTHON + ' python_scripts/concat_YZ_mpi.py '
	command = cmd_base + root_dir + ' ' + str(p1) + ' ' + str(p2) + ' ' + str(start_slice) + ' ' +  str(end_slice) + ' ' + str(inc) + ' ' + str(iproc) + ' ' + fnrs
	os.system(command)
	if( delete_XY ):
		command = 'rm -f ' + root_dir + 'output/2D/' + fnrs + '_*.nc'
		os.system(command)
	command = 'ls -lh ' + root_dir + 'output/slices/2D/' + fnrs + '* > paths/concatenated_' + fnrs + '_files'
	os.system(command)



#---------------------------------------------------- 
# Task 3: concatenate the XY plane nc files together
# creating a separate file for each time slice.
# Results stored in output/slices/2D.
# e.g.
# Saved filenames: fnrs_iproc-jproc.nc
# XY_000-001.nc	XY_001-001.nc
#----------------------------------------------------
if( do_XY ):
	fnrs = 'XY'                             # filename root string
	start_slice = 0                         # starting time slice  i.e. slices = np.arange(start_slice,end_slice,inc)
	end_slice = 51                          # ending time slice    
	inc = 1                                 # time slice increment
	jproc = 1                               # jproc for the XY plane (depends on z level of XY plane saved)
	numprocs = 4                            # number of processors to be used by concat_XY_mpi.py	
	
	cmd_base = MPIRUN + ' -np ' + str(numprocs) + ' ' + PYTHON + ' python_scripts/concat_XY_mpi.py '
	command = cmd_base + root_dir + ' ' + str(p1) + ' ' + str(p2) + ' ' + str(start_slice) + ' ' + str(end_slice) + ' ' + str(inc) + ' ' + str(iproc) + ' ' + fnrs
	os.system(command)
	if( delete_XY ):
		command = 'rm -f ' + root_dir + 'output/2D/' + fnrs + '_*.nc'
		os.system(command)
	command = 'ls -lh ' + root_dir + 'output/slices/2D/' + fnrs + '* > paths/concatenated_' + fnrs + '_files'
	os.system(command)


#---------------------------------------------------- 
# Task 4: concatenate the 3d XYZ nc files together.
# Results stored in output/slices/3D.
#----------------------------------------------------
if( do_XYZ ):
	fnrs = 'XYZ'                            # filename root string
	start_slice = 0                         # starting time slice  i.e. slices = np.arange(start_slice,end_slice,inc)
	end_slice = 1000                        # ending time slice    
	inc = 200                               # time slice increment
	numprocs = 5                            # number of processors to be used by concat_YZ_mpi.py	
	
	cmd_base = MPIRUN + ' -np ' + str(numprocs) + ' ' + PYTHON + ' python_scripts/concat_XYZ_mpi.py '
	command = cmd_base + root_dir + ' ' + str(p1) + ' ' + str(p2) + ' ' + str(start_slice) + ' ' + str(end_slice) + ' ' + str(inc) + ' ' + fnrs
	os.system(command)
	if( delete_XYZ ):
		command = 'rm -f ' + root_dir + 'output/3D/' + fnrs + '_*.nc'
		os.system(command)
	command = 'ls -lh ' + root_dir + 'output/slices/3D/' + fnrs + '* > paths/concatenated_' + fnrs + '_files'
	os.system(command)



