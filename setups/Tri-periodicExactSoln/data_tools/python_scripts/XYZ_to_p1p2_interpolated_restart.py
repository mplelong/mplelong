#!python
"""
    Script to create a set of restart files from a concatenated XYZ_xxxxxx.nc data file.
    
    NX, NY and NZ are the desired number of global gridpoints
    p1 and p2 are the decomposition parameters for the restart run, it doesn't matter
    how many were used in the original run.
    
    Both the input global 3d file and the distributed output files have the
    storage order (t,z,y,x) which matches that of the distributed 3d files
    written by flow_solve.
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,sys
import numpy as np
import scipy as sp
from scipy.io import netcdf  
#-----------------------------------------------------------------------------------------

#-------------------------
# path to python 
#-------------------------
PYTHON="$PYTHON"    # get env variable



#--------------------------------------------------------------------------- 
# top directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
top_dir = sys.argv[1]
#  add / for convenience
top_dir = top_dir + '/'
restart_dir = top_dir + 'RESTART/'           # where the p1xp2 restart files will be written, created if necessary
slice_dir = top_dir + 'output/slices/3D/'    # where the global 3d input file lives

# input file
infile = sys.argv[2]
fname = slice_dir + infile

# global number of grid points to interpolate to
NX= int(sys.argv[3])
NY= int(sys.argv[4])
NZ= int(sys.argv[5])

# decomposition parameters to be used in restarted flow_solve run
p1 = int(sys.argv[6]) 
p2 = int(sys.argv[7])

interp_method = sys.argv[8]

print('...     XYZ_to_p1p2_interpolated_restart.py executing... root_dir: ', top_dir)
print('...     do_second_scalar assumed FALSE, [s1]=[kg/m3] ')
print('...     input filename: ',fname)
print('...     interpolation method: ',interp_method)
print('...     output global resolution NX x NY x NZ: ',NX,NY,NZ)
print('...     output decomposition p1xp2: ',p1,p2)
print('...     output file directory:  ',restart_dir)
print('...     output filename root:  restart_')


# (a) interpolate in x and y
command = PYTHON + ' input/data_tools/python_scripts/interp_H.py ' + top_dir + ' ' + infile + ' ' + str(NX) + ' ' + str(NY) + ' ' + interp_method
print( command )
os.system(command)
print(' ')


# (b) interpolate in z
command = PYTHON + ' input/data_tools/python_scripts/interp_V.py ' + top_dir + ' ' + 'horiz_interp.nc XYZ_interpolated.nc' + ' ' + str(NZ) + ' ' + interp_method
print( command )
os.system(command)
print(' ')


# (c) decompose interpolated result into p1xp2 restart files
command = PYTHON + ' input/data_tools/python_scripts/XYZ_to_p1p2_restart.py ' + top_dir + ' ' + str(p1) + ' ' + str(p2) + ' XYZ_interpolated.nc'
print( command )
os.system(command)
print(' ')
print(' ')
print(' ')


print('... Cleaning up ')
command = 'rm -f output/slices/3D/horiz_interp.nc output/slices/3D/XYZ_interpolated.nc'
print( command )
os.system(command)


print('... Done. ')
print('...       The restart files are in ',restart_dir)
command = 'ls -lh ' + restart_dir
print('          ', command )
os.system(command)
