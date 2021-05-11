#!/opt/local/bin/python2.7
"""
    Script to create a set of restart files from a concatenated yz_slice_xxxxxx.nc file.
    NP is the number of processors to be used for the restarted run, it doesn't matter
    how many were used in the original run.
    The resulting files are (by convention) written to flow_solve/RESTART directory
    with root names of restart_

    MODIFIED VERSION SO THAT THE OUTPUT FILES HAVE THE SAME INDEXING CONVENTION
    AS DISTRIBUTED XYZ FILES WRITTEN BY FLOW_SOLVE: (t,z,y,x)
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os
import numpy as np
import scipy as sp
from scipy.io import netcdf                   # my usual way of reading
from netCDF4 import Dataset        # use this to avoid scipy error in writing unrecognizable netcdf format
#-----------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#  params pointing to global file w/ initial conditions data
#---------------------------------------------------------------------------------------
top_dir = '/Users/kraig/kbw_git/flow_solve/'
islice = 20    # 
do_second_scalar=False

slice_dir = top_dir + 'output/slices/2D/YZ/'
Lx = 1.0                # [m] width of domain for restart
nx = 1                  # number of x slices in global restart data file
add_x_slices = 0        # whether to distribute 2d data into 3d, i.e. in x (periodic)

p1 = 1                  # decomposition parameters for restart
p2 = 4                  # decomposition parameters for restart
if( add_x_slices==1 ):
    nx = 96             # desired number of global x pts for 3d restart from a 2d (YZ) run
    Lx = 0.25           # [m] 

dx = Lx/nx              # [m] dx for restart file
X = np.linspace(0,Lx-dx,nx)


#---------------------------------------------------------------------------------------
#  read in data from the desired z-concatenated netcdf file
#---------------------------------------------------------------------------------------
fname = slice_dir + 'YZ_' + str(islice).zfill(6) + '.nc'
ncf = netcdf.netcdf_file(fname, 'r')
Y = ncf.variables['y'].data.copy()      # [m]
Z = ncf.variables['z'].data.copy()      # [m]
U = ncf.variables['u'].data.copy()      # [m/s]
V = ncf.variables['v'].data.copy()      # [m/s]
W = ncf.variables['w'].data.copy()      # [m/s]
S1 = ncf.variables['s1'].data.copy()    # [kg/m3]
if( do_second_scalar):
    S2 = ncf.variables['s2'].data.copy()    # [1]
T = ncf.variables['time'].data.copy()   # [s]
ncf.close()
print " restart time in seconds:  ",T

ny = Y.size
nz = Z.size
print np.shape(V)

outdir = top_dir + 'RESTART'
if not os.path.exists(outdir):
    os.makedirs(outdir)

#---------------------------------------------------------------------------------------
#  number of z points per processor (sin/cos expansions), last processor gets 1 extra
#--------------------------------------------------------------------------------------
locnz = (nz-1)/p2
print "locnz for all but the last processor: %d" % (locnz)
#---------------------------------------------------------------------------------------
#  number of x points per processor (fourier expansions)
#--------------------------------------------------------------------------------------
if( add_x_slices==1 ):
    locnx = (nx)/p1
else:
    locnx = 1
print "locnx                               : %d" % (locnx)

for jproc in range(p2):
 k0 = jproc*locnz
 k1 = k0 + locnz
 if jproc==p2-1:
    locnz = (nz-1)/p2 + 1
    k1 = k1 + 1

 for iproc in range(p1):

    i0 = iproc*locnx
    i1 = i0 + locnx

    fname = outdir + '/restart_' + str(iproc).zfill(3) + '-' + str(jproc).zfill(3) + '.nc'

    #f = netcdf.netcdf_file(fname, 'w')
    f = Dataset(fname, "w", format="NETCDF4")
    print "created output netcdf file for j processor: %d" % (jproc)

    f.history = 'Created by create_restart_files.py'
    f.createDimension('idimension', locnx)
    f.createDimension('jdimension', ny)
    f.createDimension('kdimension', locnz)
    f.createDimension('timedimension', 1)
    print "-- Created dimensions"

    #----------------------------------------------------------------------
    #               name datatype dimensions      'd' is double precision
    #----------------------------------------------------------------------
    x  = f.createVariable('x', 'd', ('idimension',))
    y  = f.createVariable('y', 'd', ('jdimension',))
    z  = f.createVariable('z', 'd', ('kdimension',))
    time = f.createVariable('time', 'd', ('timedimension',))
    print "-- Created dimension variables"

    u  = f.createVariable('u', 'd', ('timedimension','kdimension','jdimension','idimension'))
    v  = f.createVariable('v', 'd', ('timedimension','kdimension','jdimension','idimension'))
    w  = f.createVariable('w', 'd', ('timedimension','kdimension','jdimension','idimension'))
    s1 = f.createVariable('s1','d', ('timedimension','kdimension','jdimension','idimension'))
    if( do_second_scalar):
        s2 = f.createVariable('s2','d', ('timedimension','kdimension','jdimension','idimension'))
    print "-- Created variables"

    time.units = 's'
    x.units = 'm'
    y.units = 'm'
    z.units = 'm'
    u.units = 'm/s'
    v.units = 'm/s'
    w.units = 'm/s'
    s1.units = 'deg C'
    if( do_second_scalar):
        s1.units = '1'
    print "-- attached units strings"

    time[:] = T
    print np.shape(x),np.shape(X[i0:i1])
    if( add_x_slices==1 ):
        x[:] = X[i0:i1]
    else:
        x[:] = X[:]

    y[:] = Y[:]
    print np.shape(z),np.shape(Z[k0:k1])
    z[:] = Z[k0:k1]

    t_id = 0   # only time slice
    ii = 0     # 1st (possibly only) x slice
    print np.shape(v),np.shape(V[k0:k1,t_id,:])
    u[0,:,:,ii]  =  U[k0:k1,t_id,:]
    v[0,:,:,ii]  =  V[k0:k1,t_id,:]
    w[0,:,:,ii]  =  W[k0:k1,t_id,:]
    s1[0,:,:,ii] = S1[k0:k1,t_id,:]
    if( do_second_scalar):
        s2[0,:,:,ii] = S2[k0:k1,t_id,:]

    if( add_x_slices == 1 ):
        for ii in range(locnx):
            u[0,:,:,ii]  =  u[0,:,:,0]
            v[0,:,:,ii]  =  v[0,:,:,0]
            w[0,:,:,ii]  =  w[0,:,:,0]
            s1[0,:,:,ii]  =  s1[0,:,:,0]
            if( do_second_scalar):
                s2[0,:,:,ii]  =  s2[0,:,:,0]
    print "-- filled/wrote data values"

    f.close()
