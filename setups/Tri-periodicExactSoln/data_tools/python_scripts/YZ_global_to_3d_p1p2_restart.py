#!python
"""
    Script to create a set of restart files from a concatenated yz_slice_xxxxxx.nc file.
    p1 and p2 are the decomposition parameters for the restarted run, it doesn't matter
    how many were used in the original run.
    
    MODIFIED VERSION SO THAT THE OUTPUT FILES HAVE THE SAME INDEXING CONVENTION
    AS DISTRIBUTED XYZ FILES WRITTEN BY FLOW_SOLVE: (t,z,y,x)
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,sys
import numpy as np
import scipy as sp
from scipy.io import netcdf                   # my usual way of reading
#-----------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#  params pointing to global file w/ initial conditions data
#---------------------------------------------------------------------------------------
top_dir = '../../../'                        # e.g. the top flow_solve directory; scripts live in input/data_tools
slice_dir = top_dir + 'output/slices/2D/'    # where the concatenated 2d slice lives
restart_dir = top_dir + 'input/RESTART/'     # where the p1xp2 restart files will be written, created if necessary
islice = 55                                  # time slice string in global YZ filename
do_second_scalar=True                        # whether or not s2 data is stored in the input/output files
do_U=True                                    # whether or not u data is stored in the input/output files
add_x_slices=False                           # whether to distribute 2d data into 3d, i.e. in x (periodic)
p1 = 1                                       # decomposition parameters for restart
p2 = 4                                       # decomposition parameters for restart


#---------------------------------------------------------------------------------------
#  set up x coordinate info
#---------------------------------------------------------------------------------------
if( add_x_slices ):
	nx = 385            # desired number of global x pts for restart
	Lx = 0.76           # [m]
else:
	nx = 1              # default for 2d in YZ plane
	Lx = 1.0            # [m]  default for 2d in YZ plane    
dx = Lx/nx              # [m] dx for restart file
X = np.linspace(0,Lx-dx,nx)


#---------------------------------------------------------------------------------------
#  read in data from the desired z-concatenated 2d (y-z) netcdf file
#---------------------------------------------------------------------------------------
fname = slice_dir + 'YZ_' + str(islice).zfill(6) + '.nc'
ncf = netcdf.netcdf_file(fname, 'r')
Y = ncf.variables['y'].data.copy()          # [m]   
Z = ncf.variables['z'].data.copy()          # [m]
if( do_U ): 
	U = ncf.variables['u'].data.copy()      # [m/s]
V = ncf.variables['v'].data.copy()          # [m/s] 
W = ncf.variables['w'].data.copy()          # [m/s]
S1 = ncf.variables['s1'].data.copy()        # [kg/m3]
if( do_second_scalar):
    S2 = ncf.variables['s2'].data.copy()    # [1]
T = ncf.variables['time'].data.copy()       # [s]
ncf.close()

ny = Y.size
nz = Z.size
print np.shape(V)    # (nz,nt=1,ny)


outdir = restart_dir
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
#---------------------------------------------------------------------------------------
#  number of z points per processor (sin/cos expansions), last processor gets 1 extra
#--------------------------------------------------------------------------------------
locnz = (nz-1)/p2
print "locnz for all but the last processor: %d" % (locnz)
#---------------------------------------------------------------------------------------
#  number of x points per processor (assuming fourier expansions)
#--------------------------------------------------------------------------------------
if( add_x_slices ):
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
    
    f = netcdf.netcdf_file(fname, 'w')
    print "created output netcdf file for j processor: %d" % (jproc)
    
    f.history = 'Created by YZ_global_to_3d_p1p2_restart.py'
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
    s1.units = 'kg/m3'
    if( do_second_scalar):
        s1.units = '1'
    print "-- attached units strings"
    
    time[:] = T
    print np.shape(x),np.shape(X[i0:i1])
    if( add_x_slices ):
        x[:] = X[i0:i1] 
    else:
        x[:] = X[:]
    
    y[:] = Y[:]
    print np.shape(z),np.shape(Z[k0:k1])
    z[:] = Z[k0:k1]
    
    t_id = 0   # only time slice
    ii = 0     # 1st (possibly only) x slice
    print np.shape(u),np.shape(U[k0:k1,t_id,:])
    u[0,:,:,ii]  =  U[k0:k1,t_id,:]
    v[0,:,:,ii]  =  V[k0:k1,t_id,:]
    w[0,:,:,ii]  =  W[k0:k1,t_id,:]
    s1[0,:,:,ii] = S1[k0:k1,t_id,:]
    if( do_second_scalar):
        s2[0,:,:,ii] = S2[k0:k1,t_id,:]
    
    if( add_x_slices ):
        for ii in range(locnx):
            u[0,:,:,ii]  =  u[0,:,:,0]
            v[0,:,:,ii]  =  v[0,:,:,0]
            w[0,:,:,ii]  =  w[0,:,:,0]
            s1[0,:,:,ii]  =  s1[0,:,:,0]
            if( do_second_scalar):
                s2[0,:,:,ii]  =  s2[0,:,:,0]
    print "-- filled/wrote data values"
    
    f.close()
