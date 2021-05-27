#!python
"""
    Script to read in a global XYZ netcdf file and interpolate in both
    x and y, writing the results to a new global XYZ netcdf file.
    Both input and output files have data layout (t,z,y,x) according to ncdump -h
    this matches the layout of the  distributed 3d files written by flow_solve.
    
    Script assumes the input file lives in $FLOW_SOLVE_ROOT/output/slices/3D.
    The output file is written to horiz_interp.nc in the same directory as the input file
    
    Modified to reflect that x is periodic, y and z are not in this problem
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,sys
import numpy as np
import scipy as sp
from scipy.io import netcdf                 
from scipy import interpolate

#--------------------------------------------------------------------------- 
# flow_solve root directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
root_dir = sys.argv[1]
#  add / for convenience
root_dir = root_dir + '/'

infile = sys.argv[2]    # input XYZ netcdf filename e.g. XYZ_000768.nc
infile_dir = root_dir + 'output/slices/3D/'

NX = int(sys.argv[3])
NY = int(sys.argv[4])
interp_method = sys.argv[5]
print('...     interp_H executing: root_dir, infile, requested NX and NY: ', root_dir,infile)
print('...     do_second_scalar assumed FALSE, [s1]=[kg/m3] infile in ',infile_dir)
print('...     all y and z  CLOSED      x PERIODIC                    ')
print('...     output dir and filename: ',root_dir,'    output/slices/3D/horiz_interp.nc ')
print('...     XY interpolation method:  ',interp_method,'   NX,NY = ',NX,NY)


#---------------------------------------------------------------------------------------
#  params pointing to global XYZ file at initial resolution
#---------------------------------------------------------------------------------------
slice_dir = root_dir + 'output/slices/3D/' 
fname1 = slice_dir + infile                   # original global XYZ file
fname2 = slice_dir + 'horiz_interp.nc'        # output global XYZ file
do_second_scalar=False                        # s2 stored in original file?

#---------------------------------------------------------------------------------------
#  read in independent variable data from input netcdf file
#---------------------------------------------------------------------------------------
ncf = netcdf.netcdf_file(fname1, 'r',version=2)
T = ncf.variables['time'].data.copy()   # [s]
Y = ncf.variables['y'].data.copy()      # [m]
Z = ncf.variables['z'].data.copy()      # [m]
X = ncf.variables['x'].data.copy()      # [m]

#---------------------------------------------------------------------------------------
#  shift grid so origin is at zero....convention for nested restart
# (true anyway for global XYZ files)
#---------------------------------------------------------------------------------------
X = X - X[0]
Y = Y - Y[0]


#---------------------------------------------------------------------------------------
#  original grid, assume *CLOSED* domain in y,z  PERIODIC in x
#---------------------------------------------------------------------------------------
dx = X[1]-X[0]       ; dy = Y[1]-Y[0]  ; dz = Z[1]-Z[0]
Lx = X[-1]-X[0]+dx   ; Ly = Y[-1]-Y[0] ; Lz = Z[-1]-Z[0]
nx = X.size ; ny = Y.size; nz = Z.size

#---------------------------------------------------------------------------------------
# grid at new resolution  X is open interval, Y is closed interval
#---------------------------------------------------------------------------------------
XI = np.linspace(0, Lx-Lx/NX, NX) 
XI = np.arange(0,NX)*(Lx/NX) 
YI = np.linspace(0, Ly, NY)



#---------------------------------------------------------------------------------------
# create the output file  (fname2)
#---------------------------------------------------------------------------------------
f = netcdf.netcdf_file(fname2, 'w',version=2)
f.createDimension('idimension', NX)
f.createDimension('jdimension', NY)
f.createDimension('kdimension', nz)
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
x[0:NX] = XI[0:NX]
y[0:NY] = YI[0:NY]
z[0:nz] = Z[0:nz]

for k in range(0, nz):
   if(np.mod(k,10)==0): print('k =', k)   
   U = ncf.variables['u'][0,k,0:ny,0:nx].copy()      # [m/s]   read kth xy plane
   V = ncf.variables['v'][0,k,0:ny,0:nx].copy()      # [m/s]
   W = ncf.variables['w'][0,k,0:ny,0:nx].copy()      # [m/s]
   S1 = ncf.variables['s1'][0,k,0:ny,0:nx].copy()    # [kg/m3]
   U = U.squeeze()
   V = V.squeeze()
   W = W.squeeze()
   S1 = S1.squeeze()
   g = interpolate.interp2d(X, Y, U,interp_method)   
   UI = g(XI, YI)   
   g = interpolate.interp2d(X, Y, V,interp_method)
   VI = g(XI, YI)
   g = interpolate.interp2d(X, Y, W,interp_method)
   WI = g(XI, YI)
   g = interpolate.interp2d(X, Y, S1,interp_method)
   S1I = g(XI, YI)
   u[0,k,0:NY,0:NX]  =  UI[0:NY,0:NX]
   v[0,k,0:NY,0:NX]  =  VI[0:NY,0:NX]
   w[0,k,0:NY,0:NX]  =  WI[0:NY,0:NX]
   s1[0,k,0:NY,0:NX] =  S1I[0:NY,0:NX]
   
   
   if( do_second_scalar ):
   		S2 = ncf.variables['s2'][0,k,0:ny,0:nx].copy()    # 
   		S2 = S2.squeeze()
   		g = interpolate.interp2d(X, Y, S2,interp_method)
   		S2I = g(XI, YI)
   		s2[0,k,0:NY,0:NX] =  S2I[0:NY,0:NX]
   		
ncf.close()
f.close()
