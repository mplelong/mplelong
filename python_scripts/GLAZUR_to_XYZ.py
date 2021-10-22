#!python
"""
    Script to create a a concatenated XYZ_xxxxxx.nc data file 
    from GLAZUR output to match the convention in Flow_Solve. 
    This creates the concatenated XYZ_xxxxxx.n  file
    storage order (t,z,y,x) matches that of the distributed 3d files
    written by flow_solve

    arg1[] is top directory
    arg2[] is name of GLAZUR64 output
    arg3[] defines the x,y and z axis names e.g. 2_258 
    
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
GLAZUR_dir = top_dir +'GLAZUR_output/' # where the GLAZUR output resides
XYZ_dir = top_dir + 'XYZ/'            # where the concatenated XYZ file will be written

# input file
infile = sys.argv[2]
fname = GLAZUR_dir + infile

print('...     GLAZUR_to_XYZ.py executing... root_dir: ', top_dir)
print('...     do_second_scalar assumed FALSE, [s1]=[kg/m3] ')
print('...     input filename: ',fname)
print('...     output file directory:  ',XYZ_dir)
print('...     output filename root:  XYZ_')


#------------------------------------------------------------------------
# path to python w/ mpi4py and corresponding mpirun
#------------------------------------------------------------------------
PYTHON="$PYTHON"    # get env variable
MPIRUN="$MPIRUN"    # get env variable



#---------------------------------------------------------------------------------------
#  params pointing to input netcdf file w/  global 3D fields at one time
#---------------------------------------------------------------------------------------

#islice = 0                                   # typically only 1 time slice per XYZ file
do_second_scalar=False                       # s2 present in original file??

#---------------------------------------------------------------------------------------
#  desired output directory for concatenated files
#---------------------------------------------------------------------------------------
outdir = XYZ_dir
if not os.path.exists(outdir):
    os.makedirs(outdir)



#---------------------------------------------------------------------------------------
#  read in the independent variable data, need to set loop limits etc
#---------------------------------------------------------------------------------------
# This is the only thing I've tried that has worked for reading in axis names
X = str('X_AXIS'+ sys.argv[3])
Y = str('Y_AXIS' + sys.argv[3])
Z = str('Z_AXIS' + sys.argv[3])

ncf = netcdf.netcdf_file(fname, 'r',version=2)
x_GL = ncf.variables[X].data[:].copy()                               # [m]
y_GL = ncf.variables[Y].data[:].copy()                               # [m]
z_GL = ncf.variables[Z].data[:].copy()                               # [m]
nx=len(x_GL); ny=len(y_GL); nz=len(z_GL)
print('size of global 3D data nx,ny,nz:   ', nx,ny,nz)
time = ncf.variables['TIME'].data[0].copy()                            # [s]
print('time stamp for global 3D data file:  ',time,'  [s]')

# read in the dependent variable data

u_GL = ncf.variables['U'].data[:].copy()
v_GL = ncf.variables['V'].data[:].copy()
w_GL = ncf.variables['W'].data[:].copy()
s1_GL = ncf.variables['S1'].data[:].copy()
rho_bar=ncf.variables['S1_BAR'].data[:].copy()
rho_bar_z=ncf.variables['S1_BAR_Z'].data[:].copy()
rho_bar_zz=ncf.variables['S1_BAR_ZZ'].data[:].copy()

print('have read in the dependent variables')
# create the output file in flow_solve format    
fname = outdir + 'XYZ_000000.nc'                
f = netcdf.netcdf_file(fname, 'w')
print (fname)
f.history = 'Created by create_restart_files.py'
f.createDimension('idimension', nx)
f.createDimension('jdimension', ny)
f.createDimension('kdimension', nz)
f.createDimension('timedimension', 1)
f.createDimension('s1bardimension',3)
print ("-- Created dimensions",nx,ny,nz)
    
        #----------------------------------------------------------------------
        #               name datatype dimensions      'd' is double precision
        #----------------------------------------------------------------------
x  = f.createVariable('x', 'd', ('idimension',))
y  = f.createVariable('y', 'd', ('jdimension',))
z  = f.createVariable('z', 'd', ('kdimension',))
time = f.createVariable('time', 'd', ('timedimension',))
print ("-- Created dimension variables")
        
u  = f.createVariable('u', 'd', ('timedimension','kdimension','jdimension','idimension'))
v  = f.createVariable('v', 'd', ('timedimension','kdimension','jdimension','idimension'))
w  = f.createVariable('w', 'd', ('timedimension','kdimension','jdimension','idimension'))
s1 = f.createVariable('s1','d', ('timedimension','kdimension','jdimension','idimension'))
print("-- Created u,v,w and s1 variables")

s1_bar = f.createVariable('s1_bar','d',('kdimension','s1bardimension'))
print("-- Created s1_bar variable")

if( do_second_scalar):
    s2 = f.createVariable('s2','d', ('timedimension','kdimension','jdimension','idimension'))

print ("-- Finished with all variable creation")
    
time.units = 's'
x.units = 'm'
y.units = 'm'
z.units = 'm'
u.units = 'm/s'
v.units = 'm/s'
w.units = 'm/s'
s1.units = 'kg/m3'

# not sure how to attribute dimension to s1_bar array

if( do_second_scalar):
   s1.units = '1'
print ("-- attached units strings")

print(np.shape(time))    
time[:] = time[:] 
print (np.shape(x))
x[:] = x_GL[:] 
    
print (np.shape(y),np.shape(y))    
y[:] = y_GL[:]
print (np.shape(z),np.shape(z))
z[:] = z_GL[:]   
print (np.shape(u))

        
        #----------------------------------------------------------------------
        #  write data into fname with the time slice equal to zero
        #----------------------------------------------------------------------    
u[0,:,:,:]  =  u_GL[:]
v[0,:,:,:]  =  v_GL[:]
w[0,:,:,:]  =  w_GL[:]
s1[0,:,:,:] =  s1_GL[:]
s1_bar[:,0] = rho_bar[:]
s1_bar[:,1] = rho_bar_z[:]
s1_bar[:,2] = rho_bar_zz[:]
if( do_second_scalar):
   s2[0,:,:,:] = s2
 
f.close() ; del u,v,w,s1
if( do_second_scalar): del s2

