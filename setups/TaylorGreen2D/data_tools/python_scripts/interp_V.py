#!python
"""
    Script to read in a global XYZ netcdf file and interpolate in z,
    writing the results to a new global XYZ netcdf file.
    Both input and output files have data layout (t,z,y,x) according to ncdump -h
    This matches the layout of the  distributed 3d files written by flow_solve. 
        
    Script assumes the input and output files are in $FLOW_SOLVE_ROOT/output/slices/3D
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,sys
import numpy as np
import scipy as sp
from scipy.io import netcdf                   # my usual way of reading
from scipy import interpolate

#--------------------------------------------------------------------------- 
# flow_solve root directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
root_dir = sys.argv[1]
#  add / for convenience
root_dir = root_dir + '/'

infile  = sys.argv[2]    # input XYZ netcdf filename 
outfile = sys.argv[3]    # output XYZ netcdf filename 
infile_dir = root_dir + 'output/slices/3D/'

NZ = int(sys.argv[4])
interp_method = sys.argv[5]
print('...     interp_V executing... root_dir: ', root_dir)
print('...     do_second_scalar assumed FALSE, [s1]=[kg/m3] ')
print('...     z coordinate dimension interval assumed to be CLOSED ')
print('...     input and output directory: ',infile_dir)
print('...     input and output filenames: ',infile,outfile)
print('...     Z interpolation method:  ',interp_method,'   NZ = ',NZ)

#---------------------------------------------------------------------------------------
#  params pointing to global XYZ file 
#---------------------------------------------------------------------------------------
slice_dir = root_dir + 'output/slices/3D/'  
fname1 = slice_dir + infile         # input XYZ file
fname2 = slice_dir + outfile        # output XYZ file
do_second_scalar=False              # s2 stored in original file?

#---------------------------------------------------------------------------------------
#  read in data from the input file  (fname1)
#---------------------------------------------------------------------------------------
print("opening input data file: ",fname1)
ncf = netcdf.netcdf_file(fname1, 'r',version=2)
T = ncf.variables['time'].data.copy()   # [s]
Y = ncf.variables['y'].data.copy()      # [m]
Z = ncf.variables['z'].data.copy()      # [m]
X = ncf.variables['x'].data.copy()      # [m]

#---------------------------------------------------------------------------------------
#  shift grid so origin is at zero....convention for nested restart
# (true anyway for global XYZ files)
#---------------------------------------------------------------------------------------
Z = Z - Z[0]    # X Y already shifted in interp_H.py

#---------------------------------------------------------------------------------------
#  original grid, assume *CLOSED* domain in z  
#---------------------------------------------------------------------------------------
dz = Z[1]-Z[0]
Lz = Z[-1]-Z[0]
nx = X.size ; ny = Y.size; nz = Z.size


#---------------------------------------------------------------------------------------
# grid at new resolution 
#---------------------------------------------------------------------------------------
ZI = np.linspace(0, Lz, NZ)



#---------------------------------------------------------------------------------------
# create the output file  (fname2)
#---------------------------------------------------------------------------------------
print("opening output data file: ",fname2)
f = netcdf.netcdf_file(fname2,'w',version=2)
f.createDimension('idimension', nx)
f.createDimension('jdimension', ny)
f.createDimension('kdimension', NZ)
f.createDimension('timedimension', 1)
#----------------------------------------------------------------------
#               name datatype dimensions      'd' is double precision
#----------------------------------------------------------------------
x  = f.createVariable('x', 'd', ('idimension',))
y  = f.createVariable('y', 'd', ('jdimension',))
z  = f.createVariable('z', 'd', ('kdimension',))
time = f.createVariable('time', 'd', ('timedimension',))
u  = f.createVariable('u', 'd', ('timedimension', 'kdimension', 'jdimension', 'idimension'))
v  = f.createVariable('v', 'd', ('timedimension', 'kdimension', 'jdimension', 'idimension'))
w  = f.createVariable('w', 'd', ('timedimension', 'kdimension', 'jdimension', 'idimension'))
s1 = f.createVariable('s1','d', ('timedimension', 'kdimension', 'jdimension', 'idimension'))
if( do_second_scalar ):
	s2 = f.createVariable('s2','d', ('timedimension', 'kdimension', 'jdimension', 'idimension'))
	s2.units = 'dless'
time.units = 's'
x.units = 'm'
y.units = 'm'
z.units = 'm'
u.units = 'm/s'
v.units = 'm/s'
w.units = 'm/s'
s1.units = 'kg/m3'
time[:] = T
x[:] = X[:]
y[:] = Y[:]
z[0:NZ] = ZI[0:NZ]

for j in range(0, ny):
   if(np.mod(ny,10)==0): print("j = ", j)
   for i in range(0, nx): 
      U = ncf.variables['u'][0,0:nz,j,i].copy().squeeze()      # [m/s]
      V = ncf.variables['v'][0,0:nz,j,i].copy().squeeze()      # [m/s]
      W = ncf.variables['w'][0,0:nz,j,i].copy().squeeze()      # [m/s]
      S1 = ncf.variables['s1'][0,0:nz,j,i].copy().squeeze()    # [kg/m3] 
       
      #   cubic throws an MKL library error on Atlantic        
      g = interpolate.interp1d(Z, U, kind=interp_method)
      UI = g(ZI)
      g = interpolate.interp1d(Z, V, kind=interp_method)
      VI = g(ZI)
      g = interpolate.interp1d(Z, W, kind=interp_method)
      WI = g(ZI)
      g = interpolate.interp1d(Z, S1,kind=interp_method)
      S1I = g(ZI)
      
      u[0,0:NZ,j,i]  =  UI[0:NZ]
      v[0,0:NZ,j,i]  =  VI[0:NZ]
      w[0,0:NZ,j,i]  =  WI[0:NZ]
      s1[0,0:NZ,j,i] =  S1I[0:NZ]
      
      if( do_second_scalar ):
      	S2 = ncf.variables['s2'][0,0:nz,j,i].copy()
      	S2 = S2.squeeze()
      	g = interpolate.interp1d(Z, S2, kind=interp_method)
      	S2I = g(ZI)
      	s2[0,0:NZ,j,i] =  S2I[0:NZ]
      	

ncf.close()
f.close()
