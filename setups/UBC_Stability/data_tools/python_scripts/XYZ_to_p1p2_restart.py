#!python
"""
    Script to create a set of restart files from a concatenated XYZ_xxxxxx.nc data file.
    
    p1 and p2 are the decomposition parameters for the restarted run, it doesn't matter
    how many were used in the original run.
    
    both the input global 3d file and the distributed output files have the
    storage order (t,z,y,x) which matches that of the distributed 3d files
    written by flow_solve
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
restart_dir = top_dir + 'RESTART/'           # where the p1xp2 restart files will be written, created if necessary
slice_dir = top_dir + 'output/slices/3D/'    # where the global 3d file lives

# decomposition parameters to be used in restarted flow_solve run
p1 = flow_solve_root = int(sys.argv[2]) 
p2 = flow_solve_root = int(sys.argv[3])

# input file
infile = sys.argv[4]
fname = slice_dir + infile

print('...     XYZ_to_p1p2_restart.py executing... root_dir: ', top_dir)
print('...     do_second_scalar assumed FALSE, [s1]=[kg/m3] ')
print('...     input filename: ',fname)
print('...     output decomposition p1xp2: ',p1,p2)
print('...     output file directory:  ',restart_dir)
print('...     output filename root:  restart_')


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
X = ncf.variables['x'].data[:].copy()                               # [m]
Y = ncf.variables['y'].data[:].copy()                               # [m]
Z = ncf.variables['z'].data[:].copy()                               # [m]
nx=len(X); ny=len(Y); nz=len(Z)
print('size of global 3D data nx,ny,nz:   ', nx,ny,nz)
T = ncf.variables['time'].data[0].copy()                            # [s]
print('time stamp for global 3D data file:  ',T,'  [s]')

#---------------------------------------------------------------------------------------
#  number of z points per processor (sin/cos expansions), last processor gets 1 extra
#--------------------------------------------------------------------------------------
locnz = int((nz-1)/p2)
print("Assuming sin/cos expansions, locnz for all but the last processor: ", locnz) 

#---------------------------------------------------------------------------------------
#  number of x points per processor (cos expansions)
#--------------------------------------------------------------------------------------
locnx = int((nx-1)/p1) 
print("Assuming cos in x, locnx for but the last all processors                  : ",locnx) 


for jproc in range(p2):
    locnz = int((nz-1)/p2)
    k0 = jproc*locnz
    k1 = k0 + locnz
    if jproc==p2-1:
        locnz = int((nz-1)/p2 + 1)
        k1 = k1 + 1
    
    for iproc in range(p1):
        locnx = int((nx-1)/p1)
        i0 = iproc*locnx
        i1 = i0 + locnx
        if iproc==p1-1:
        	locnx = int((nx-1)/p1 + 1)
        	i1 = i1 + 1
    
        fname = outdir + 'restart_' + str(iproc).zfill(3) + '-' + str(jproc).zfill(3) + '.nc'               
        f = netcdf.netcdf_file(fname, 'w')
        print (fname)
        print ("created output netcdf file for j processor: ", jproc)
            
        f.history = 'Created by create_restart_files.py'
        f.createDimension('idimension', locnx)
        f.createDimension('jdimension', ny)
        f.createDimension('kdimension', locnz)
        f.createDimension('timedimension', 1)
        print ("-- Created dimensions",locnx,ny,locnz)
    
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
        if( do_second_scalar):
            s2 = f.createVariable('s2','d', ('timedimension','kdimension','jdimension','idimension'))
        print ("-- Created variables")
    
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
        print ("-- attached units strings")
    
        time[:] = T
        print (np.shape(x),np.shape(X[i0:i1]))
        x[:] = X[i0:i1] 
    
    
        y[:] = Y[:]
        print (np.shape(z),np.shape(Z[k0:k1]))
        z[:] = Z[k0:k1]   
        print (np.shape(u))

        #----------------------------------------------------------------------
        #  read the portion of data from the global file needed for this file
        #     here the XYZ global file is written (t,z,y,x)
        #----------------------------------------------------------------------
        U = ncf.variables['u'].data[0,k0:k1,:,i0:i1].copy()          # [m/s]         
        V = ncf.variables['v'].data[0,k0:k1,:,i0:i1].copy()          # [m/s]
        W = ncf.variables['w'].data[0,k0:k1,:,i0:i1].copy()          # [m/s] 
        S1= ncf.variables['s1'].data[0,k0:k1,:,i0:i1].copy()         # [kg/m3]
        if( do_second_scalar):
            S2= ncf.variables['s2'].data[0,k0:k1,:,i0:i1].copy()     # [1]
        
        
        
        #----------------------------------------------------------------------
        #  write data into fname with the time slice equal to zero
        #----------------------------------------------------------------------    
        u[0,:,:,:]  =  U
        v[0,:,:,:]  =  V
        w[0,:,:,:]  =  W
        s1[0,:,:,:] =  S1
        if( do_second_scalar):
            s2[0,:,:,:] = S2
    
        print("-- filled/wrote data values for iproc,jproc,i:",iproc,jproc,np.shape(u))
        f.close() ; del u,v,w,s1
        if( do_second_scalar): del s2

