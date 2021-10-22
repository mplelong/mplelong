#!python
"""
    Script to mean scalar profiles read from a concatenated XYZ_xxxxxx.nc data file into a separate file
    Arguments needed are: (1) top directory name and (2) name of file in top_dir+/XYZ/ directory
    
    This file will be read by flow_solve at the beginning of restart runs  
    s1_bar and (if defined) s2_bar arrays contain s1_bar, s1_bar_z, s1_bar_zz (same for s2_bar)
    
    P. Lelong (October 2021) s1_bar array in all restart files and, 
    if needed also s2_bar
    Also included is the z coordinate
. 
    This is needed if s1_bar(1-3)/s2_bar(1-3) are read from a file rather than defined in functional form. 
    
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,sys
import numpy as np
import scipy as sp
from scipy.io import netcdf  
#-----------------------------------------------------------------------------------------

#--------------------------------------------------------------------------- 
# top directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
top_dir = sys.argv[1]
#  add / for convenience
top_dir = top_dir + '/'
restart_dir = top_dir + 'RESTART/'           # where the files will be written, created if necessary
slice_dir = top_dir + 'XYZ/'    # where the global 3d file lives

# input file
infile = sys.argv[2]
fname = slice_dir + infile

print('...     s_bar.py executing... root_dir: ', top_dir)
print('...     do_second_scalar assumed FALSE, [s1]=[kg/m3] ')
print('...     input filename: ',fname)
print('...     output file directory:  ',restart_dir)


#------------------------------------------------------------------------
# path to python w/ mpi4py and corresponding mpirun
#------------------------------------------------------------------------
PYTHON="$PYTHON"    # get env variable
MPIRUN="$MPIRUN"    # get env variable



#---------------------------------------------------------------------------------------
#  params pointing to input netcdf file w/  global 3D fields at one time
#---------------------------------------------------------------------------------------

islice = 0                                   # typically only 1 time slice per XYZ file
do_second_scalar=False                       # s2 present in original file??

#---------------------------------------------------------------------------------------
#  desired output directory for restart files
#---------------------------------------------------------------------------------------
outdir = restart_dir
if not os.path.exists(outdir):
    os.makedirs(outdir)



#---------------------------------------------------------------------------------------
#  read in the independent variable data, need to set loop limits etc
#---------------------------------------------------------------------------------------
ncf = netcdf.netcdf_file(fname, 'r',version=2)
Z = ncf.variables['z'].data[:].copy()                               # [m]
nz=len(Z)
print('size of nz:   ', nz)
T = ncf.variables['time'].data[0].copy()                            # [s]
print('time stamp for global 3D data file:  ',T,'  [s]')

# read s1_bar
S1_BAR = ncf.variables['s1_bar'].data[:].copy()
n=S1_BAR.shape
print('dimensions of S1_BAR',n)
ndim=n[1]
nzglobal=n[0]
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
    
fname = outdir + 's1_bar.nc'                
f = netcdf.netcdf_file(fname, 'w')
print (fname)
print ("created output netcdf file for s1_bar")
           
f.history = 'Created by store_scalar_bar.py'
f.createDimension('ndim',3)
f.createDimension('nzglobal',nz)
print ("-- Created dimensions",ndim,nzglobal)
    
        #----------------------------------------------------------------------
        #               name datatype dimensions      'd' is double precision
        #----------------------------------------------------------------------
z  = f.createVariable('z', 'd', ('nzglobal',))
      
print ("-- Created dimension variables")

s1_bar = f.createVariable('s1_bar','d',('nzglobal','ndim'))


print ("-- Created variables")
         
z.units = 'm'
s1_bar.units = 'kg/m3,kg/m4,kg/m5'
print ("-- attached units strings")

        #----------------------------------------------------------------------
        #  write data into fname 
        #----------------------------------------------------------------------    
z[:] = Z[:]
s1_bar[:,:] = S1_BAR[:,:]
print("-- filled/wrote data values for s1_bar and z:",s1_bar,z)
f.close(); del s1_bar
if( do_second_scalar):
            S2_BAR = ncf.variables['s2_bar'].data[:].copy()
            fname2 = outdir + 's2_bar.nc'
            f2 = netcdf.netcdf_file(fname2, 'w')
            print (fname2)
            print ("created output netcdf file for s2_bar")
            f2.history = 'Created by store_scalar_bar.py'
            f2.createDimension('ndim',3)
            f2.createDimension('nzglobal',nz)
            print ("-- Created dimensions",ndim,nzglobal)
            s2_bar=f2.createVariable('s2_bar','d',('nzglobal','ndim'))
            z = f2.createVariable('z','d',('nzglobal',))
            s2_bar.units = '1'
            print ("-- attached units strings")
            z[:] = Z[:]
            s2_bar[:,:] = S2_BAR[:,:]
            print("-- filled/wrote data values for s2_bar and z:",s2_bar,z)
            f2.close() ; del s2_bar

