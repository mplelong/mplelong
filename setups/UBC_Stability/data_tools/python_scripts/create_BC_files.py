#!
"""
Routine to 
 read concatenated coarse BC data, interpolate and write a set  
 of distributed data files for a new p1xp2 fine resolution run
 
"""
#---------------------------------------------------------------------------------------
#  import and name the various python modules I'll use
#      u_dot_grad_2d(u,v,phi,Lx,Ly,xflag,yflag):
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
#from scipy.io import netcdf
import netCDF4 as nc
from mpi4py import MPI
from scipy import interpolate
#-----------------------------------------------------------------------------------------

#-------------------------
# path to python, mpirun 
#-------------------------
PYTHON="$PYTHON"    # get env variable
MPIRUN="$MPIRUN"    # get env variable

fmt = "NETCDF4"
dtype = np.float64


#--------------------------------------------------------------------------- 
# top directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
root_dir = sys.argv[1]
#  add / for convenience
root_dir = root_dir + '/'
coarse_slice_dir = root_dir + 'BVALS/'  # location of coarse, spatially global BC data files east.nc  etc
BC_file_dir = coarse_slice_dir          # location for new, fine spatial resolution, spatially distributed files

# global number of grid points for the run reading boundary files created
nx = int(sys.argv[2])
ny = int(sys.argv[3])
nz = int(sys.argv[4])

# decomposition parameters for the run reading boundary files created
# p1 is number of processors splitting nx
# p2 is number of processors splitting ny and nz
p1 = int(sys.argv[5]) 
p2 = int(sys.argv[6])

interp_method = 'cubic'
do_second_scalar = False

print('...     create_BC_files.py executing... root_dir: ', root_dir)
print('...     do_second_scalar assumed FALSE, [s1]=[kg/m3] ')
print('...     input coarse B data directory: ',coarse_slice_dir)
print('...     output global resolution nx x ny x nz: ',nx,ny,nz)
print('...     output decomposition p1xp2: ',p1,p2)
print('...     output file directory:  ',BC_file_dir)



#-----------------------------------------------------------------------------------------
#  local array sizes
#-----------------------------------------------------------------------------------------
locnx = (nx-1)/p1         # for all iprocs except p1-1 which has 1 extra
locnz = (nz-1)/p2         # for all jprocs except p2-1 which has 1 extra

#-----------------------------------------------------------------------------------------
#  calculate the global gridpoint indices in fine grid for each processor
#-----------------------------------------------------------------------------------------
i0=np.zeros(p1, dtype=int) ; i1=np.zeros(p1, dtype=int); nxvals=np.zeros(p1, dtype=int)
for iproc in range(p1):
	i0[iproc] = iproc*locnx ; i1[iproc] = i0[iproc] + locnx ; nxvals[iproc] = locnx
	if iproc==p1-1: i1[iproc] = i1[iproc] + 1 ; nxvals[iproc] = locnx + 1	

k0=np.zeros(p2, dtype=int) ; k1=np.zeros(p2, dtype=int); nzvals=np.zeros(p2, dtype=int)			
for jproc in range(p2):				
	k0[jproc] = jproc*locnz ; k1[jproc] = k0[jproc] + locnz ; nzvals[jproc] = locnz
	if jproc==p2-1: k1[jproc] = k1[jproc] + 1 ; nzvals[jproc] = locnz + 1
	
if not os.path.exists(BC_file_dir):
	cmd = 'mkdir -p ' + BC_file_dir
	os.system(cmd)


#-----------------------------------------------------------------------------------------
#  XY planes for bottom and top   after reading, storage format is u[t,y,x]
#-----------------------------------------------------------------------------------------
for loc in [0,1]:
	if(loc==0): root = 'bottom'
	if(loc==1): root = 'top'
	fname = coarse_slice_dir + root + '.nc'
	print ('Reading file ', fname)

	#---------------------------------------------------------------------------------------
	#  read in data from the desired netcdf file:  _c  means coarse quantity
	#---------------------------------------------------------------------------------------
	ncf = nc.Dataset(fname, 'r', format=fmt)
	x_c = ncf['x'][:]
	y_c = ncf['y'][:]
	t    = ncf['t'][:]
	
	# read boundary vals
	u_c  = ncf['u'][:,:,:] 
	v_c  = ncf['v'][:,:,:]
	w_c  = ncf['w'][:,:,:]
	s1_c = ncf['s1'][:,:,:] 
	if( do_second_scalar ):
		dye_c =  ncf['s2'][:,:,:]
	
	# read boundary derivs	
	u_c_z  = ncf['u_z'][:,:,:] 
	v_c_z  = ncf['v_z'][:,:,:]
	w_c_z  = ncf['w_z'][:,:,:]
	s1_c_z = ncf['s1_z'][:,:,:] 
	if( do_second_scalar ):
		dye_c_z =  ncf['s2_z'][:,:,:] 
	        
	Ly = y_c[-1]-y_c[0]                         # [m]
	Lx = x_c[-1]-x_c[0]                         # [m]
	nt = t.size
	ny_c = y_c.size
	nx_c = x_c.size
	#print( 'shape of 2d data arrays u_c, etc: ', v_c.shape, nt, ny_c, nx_c)  #[nt,ny_c,nx_c]        
	print( "...   closing ", fname)
	ncf.close()
	print (' ')
	print (' ')
	
	#---------------------------------------------------------------------------------------
	#  shift origin of coarse grid and compute fine grid mesh spacing
	#---------------------------------------------------------------------------------------
	y_c = y_c - y_c[0]
	x_c = x_c - x_c[0]
	
	dy = Ly/(ny-1.)    # dy in fine grid
	dx = Lx/(nx-1.)    # dx in fine grid
	
		
	#----------------------------------------------------------------
	#  loop over the iproc and jproc processors that will write data 
	#----------------------------------------------------------------
	if( loc==0 ): jproc = 0
	if( loc==1 ): jproc = p2-1
						
	for iproc in range(p1):
	
		y0 = 0.0 ; y1 = Ly     # y is global in YBLOCK format
		yy = np.linspace(y0,y1,ny,endpoint=True)
	
		x0 = dx*i0[iproc] ; x1 = x0 + dx*(nxvals[iproc]-1)
		xx = np.linspace(x0,x1,nxvals[iproc],endpoint=True)
		
		# vals
		U = np.zeros((nt,ny,nxvals[iproc]))
		V = np.zeros((nt,ny,nxvals[iproc]))
		W = np.zeros((nt,ny,nxvals[iproc]))
		S1 = np.zeros((nt,ny,nxvals[iproc]))
		
		# derivs
		U_z = np.zeros((nt,ny,nxvals[iproc]))
		V_z = np.zeros((nt,ny,nxvals[iproc]))
		W_z = np.zeros((nt,ny,nxvals[iproc]))
		S1_z = np.zeros((nt,ny,nxvals[iproc]))
		
		#---------------------------------------------------
		#  loop through all time slices for this processor 
		#---------------------------------------------------
		for islice in np.arange(0,nt,1):    
			
			#-------------------------------------------
			#  compute interpolating functions over
			#  entire global coarse grid 
			#-------------------------------------------
			fu  = interpolate.interp2d(x_c, y_c,  u_c[islice,:,:].squeeze(), kind='cubic')
			fv  = interpolate.interp2d(x_c, y_c,  v_c[islice,:,:].squeeze(), kind='cubic')
			fw  = interpolate.interp2d(x_c, y_c,  w_c[islice,:,:].squeeze(), kind='cubic')
			fs1 = interpolate.interp2d(x_c, y_c, s1_c[islice,:,:].squeeze(), kind='cubic')
			
			fu_z  = interpolate.interp2d(x_c, y_c,  u_c_z[islice,:,:].squeeze(), kind='cubic')
			fv_z  = interpolate.interp2d(x_c, y_c,  v_c_z[islice,:,:].squeeze(), kind='cubic')
			fw_z  = interpolate.interp2d(x_c, y_c,  w_c_z[islice,:,:].squeeze(), kind='cubic')
			fs1_z = interpolate.interp2d(x_c, y_c, s1_c_z[islice,:,:].squeeze(), kind='cubic')
		
			#-------------------------------------------------
			# get interpolated field variables
			#-------------------------------------------------			
			U[islice,:,:] = fu(xx,yy)  ; #print 'U shape ',U.shape  # ==> nt,ny,nx
			V[islice,:,:] = fv(xx,yy) 
			W[islice,:,:] = fw(xx,yy)
			S1[islice,:,:] = fs1(xx,yy)
			
			U_z[islice,:,:] = fu_z(xx,yy)  ; #print 'U shape ',U.shape  # ==> nt,ny,nx
			V_z[islice,:,:] = fv_z(xx,yy) 
			W_z[islice,:,:] = fw_z(xx,yy)
			S1_z[islice,:,:] = fs1_z(xx,yy)
			
		#-------------------------------------------------
		#  write spatially local data at all times to file
		#-------------------------------------------------				
		fname = BC_file_dir + root + '_' + str(iproc).zfill(3) + '-' + str(jproc).zfill(3) + '.nc'
		f = nc.Dataset(fname, 'w', format=fmt)
		print ('created new netcdf file: ',fname)
		timedimension = f.createDimension("timedimension", nt)
		idimension = f.createDimension("idimension", nxvals[iproc])
		jdimension = f.createDimension("jdimension", ny)
		time = f.createVariable('time', dtype, ("timedimension",))
		x  = f.createVariable('x', dtype, ("idimension",))
		y  = f.createVariable('y', dtype, ("jdimension",))
		u  = f.createVariable('u', dtype, ('timedimension','jdimension','idimension',))
		v  = f.createVariable('v', dtype, ('timedimension','jdimension','idimension',))
		w  = f.createVariable('w', dtype, ('timedimension','jdimension','idimension',))
		s1 = f.createVariable('s1',dtype, ('timedimension','jdimension','idimension',))
		if( do_second_scalar):
			s2 = f.createVariable('s2_z',dtype, ('timedimension','jdimension','idimension',))
		u_z  = f.createVariable('u_z', dtype, ('timedimension','jdimension','idimension',))
		v_z  = f.createVariable('v_z', dtype, ('timedimension','jdimension','idimension',))
		w_z  = f.createVariable('w_z', dtype, ('timedimension','jdimension','idimension',))
		s1_z = f.createVariable('s1_z',dtype, ('timedimension','jdimension','idimension',))
		if( do_second_scalar):
			s2_z = f.createVariable('s2_z',dtype, ('timedimension','jdimension','idimension',))
		
		f.history = 'Created by create_BC_files.py'				
		time.units = 's'
		y.units = 'm'
		x.units = 'm'
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		if( do_second_scalar):
			s2.units = '1'
		u_z.units = '1/s'
		v_z.units = '1/s'
		w_z.units = '1/s'
		s1_z.units = 'kg/m4'
		if( do_second_scalar):
			s2_z.units = '1/m'
				
		# write all the data to this file
		time[:] = t[:]
		x[:] =   xx[:]
		y[:] =   yy[:]
		u[:,:,:] =    U[:,:,:]
		v[:,:,:] =    V[:,:,:]
		w[:,:,:] =    W[:,:,:]
		s1[:,:,:] =  S1[:,:,:]
		u_z[:,:,:] = U_z[:,:,:]
		v_z[:,:,:] = V_z[:,:,:]
		w_z[:,:,:] = W_z[:,:,:]
		s1_z[:,:,:] = S1_z[:,:,:]
		f.close()


#-----------------------------------------------------------------------------------------
#  XZ planes for south and north   after reading, storage format is u[t,z,x]
#-----------------------------------------------------------------------------------------
for loc in [0,1]:
	if(loc==0): root = 'south'
	if(loc==1): root = 'north'
	fname = coarse_slice_dir + root + '.nc'
	print ('Reading file ', fname)

	#---------------------------------------------------------------------------------------
	#  read in data from the desired netcdf file:  _c  means coarse quantity
	#---------------------------------------------------------------------------------------
	ncf = nc.Dataset(fname, 'r', format=fmt)
	x_c = ncf['x'][:]
	z_c = ncf['z'][:]	
	t    = ncf['t'][:]
	
	# read boundary vals
	u_c  = ncf['u'][:,:,:] 
	v_c  = ncf['v'][:,:,:]
	w_c  = ncf['w'][:,:,:]
	s1_c = ncf['s1'][:,:,:] 
	if( do_second_scalar ):
		dye_c =  ncf['s2'][:,:,:]
	
	# read boundary derivs	
	u_c_y  = ncf['u_y'][:,:,:] 
	v_c_y  = ncf['v_y'][:,:,:]
	w_c_y  = ncf['w_y'][:,:,:]
	s1_c_y = ncf['s1_y'][:,:,:] 
	if( do_second_scalar ):
		dye_c_y =  ncf['s2_y'][:,:,:]    
        
	Lx = x_c[-1]-x_c[0]                         # [m]
	Lz = z_c[-1]-z_c[0]                         # [m]
	nt = t.size
	nz_c = z_c.size
	nx_c = x_c.size
	#print 'shape of 2d data arrays u_c, etc: ', v_c.shape, nt, nz_c, nx_c  #[nt,nz_c,nx_c]        
	#print "closing ", fname
	ncf.close()
	print (' ') 
	print (' ') 
	
	#---------------------------------------------------------------------------------------
	#  shift origin of coarse grid and compute fine grid mesh spacing
	#---------------------------------------------------------------------------------------
	x_c = x_c - x_c[0]
	z_c = z_c - z_c[0]
	
	dx = Lx/(nx-1.)    # dx in fine grid
	dz = Lz/(nz-1.)    # dz in fine grid
	
		
	#-------------------------------------------
	#  loop 2d iproc and jproc processors 
	#-------------------------------------------
	for iproc in range(p1):				
		for jproc in range(p2):
			
			#-------------------------------------------------
			#  determine local fine grid coordinate vals
			#-------------------------------------------------
			x0 = dx*i0[iproc] ; x1 = x0 + dx*(nxvals[iproc]-1)
			xx = np.linspace(x0,x1,nxvals[iproc],endpoint=True)
			#print x0,x1
	
			z0 = dz*k0[jproc] ; z1 = z0 + dz*(nzvals[jproc]-1)
			zz = np.linspace(z0,z1,nzvals[jproc],endpoint=True)
			#print z0,z1
			
			U = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			V = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			W = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			S1 = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			
			U_y = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			V_y = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			W_y = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			S1_y = np.zeros((nt,nzvals[jproc],nxvals[iproc]))
			
			#-------------------------------------------
			#  loop through time slices 
			#-------------------------------------------
			for islice in np.arange(0,nt,1):     # np.arange(0,nt,1):
			
				#-------------------------------------------
				#  compute interpolating functions over
				#  entire global coarse grid 
				#-------------------------------------------
				fu  = interpolate.interp2d(x_c, z_c,  u_c[islice,:,:].squeeze(), kind='cubic')
				fv  = interpolate.interp2d(x_c, z_c,  v_c[islice,:,:].squeeze(), kind='cubic')
				fw  = interpolate.interp2d(x_c, z_c,  w_c[islice,:,:].squeeze(), kind='cubic')
				fs1 = interpolate.interp2d(x_c, z_c, s1_c[islice,:,:].squeeze(), kind='cubic')
				
				fu_y  = interpolate.interp2d(x_c, z_c,  u_c_y[islice,:,:].squeeze(), kind='cubic')
				fv_y  = interpolate.interp2d(x_c, z_c,  v_c_y[islice,:,:].squeeze(), kind='cubic')
				fw_y  = interpolate.interp2d(x_c, z_c,  w_c_y[islice,:,:].squeeze(), kind='cubic')
				fs1_y = interpolate.interp2d(x_c, z_c, s1_c_y[islice,:,:].squeeze(), kind='cubic')
				
				#-------------------------------------------------
				# get interpolated field variables
				#-------------------------------------------------
				U[islice,:,:] = fu(xx,zz) ; #print 'U shape ',U.shape  # ==> nt,nz,nx     
				V[islice,:,:] = fv(xx,zz) 		
				W[islice,:,:] = fw(xx,zz)		
				S1[islice,:,:] = fs1(xx,zz)
				
				U_y[islice,:,:] = fu_y(xx,zz) ; #print 'U shape ',U.shape  # ==> nt,nz,nx     
				V_y[islice,:,:] = fv_y(xx,zz) 		
				W_y[islice,:,:] = fw_y(xx,zz)		
				S1_y[islice,:,:] = fs1_y(xx,zz)
				
			#-------------------------------------------------
			#  write local data to file
			#-------------------------------------------------				
			fname = BC_file_dir + root + '_' + str(iproc).zfill(3) + '-' + str(jproc).zfill(3) + '.nc'
			f = nc.Dataset(fname, 'w', format=fmt)
			print ('created new netcdf file: ',fname)
			timedimension = f.createDimension("timedimension", nt)
			idimension = f.createDimension("idimension", nxvals[iproc])
			kdimension = f.createDimension("kdimension", nzvals[jproc])
			time = f.createVariable('time', dtype, ("timedimension",))			
			x  = f.createVariable('x', dtype, ("idimension",))
			z  = f.createVariable('z', dtype, ("kdimension",))
			u  = f.createVariable('u', dtype, ('timedimension','kdimension','idimension',))
			v  = f.createVariable('v', dtype, ('timedimension','kdimension','idimension',))
			w  = f.createVariable('w', dtype, ('timedimension','kdimension','idimension',))
			s1 = f.createVariable('s1',dtype, ('timedimension','kdimension','idimension',))
			if( do_second_scalar):
				s2 = f.createVariable('s2_y',dtype, ('timedimension','kdimension','idimension',))
			u_y  = f.createVariable('u_y', dtype, ('timedimension','kdimension','idimension',))
			v_y  = f.createVariable('v_y', dtype, ('timedimension','kdimension','idimension',))
			w_y  = f.createVariable('w_y', dtype, ('timedimension','kdimension','idimension',))
			s1_y = f.createVariable('s1_y',dtype, ('timedimension','kdimension','idimension',))
			if( do_second_scalar):
				s2_y = f.createVariable('s2_y',dtype, ('timedimension','kdimension','idimension',))
				
			time.units = 's'
			x.units = 'm'
			z.units = 'm'
			u.units = 'm/s'
			v.units = 'm/s'
			w.units = 'm/s'
			s1.units = 'kg/m3'
			if( do_second_scalar):
				s2.units = '1'
			u_y.units = '1/s'
			v_y.units = '1/s'
			w_y.units = '1/s'
			s1_y.units = 'kg/m4'
			if( do_second_scalar):
				s2_y.units = '1/m'
					
			# write all the data to this file
			time[:] = t[:]
			z[:] =   zz[:]
			x[:] =   xx[:]
			u[:,:,:]  =    U[:,:,:]
			v[:,:,:]  =    V[:,:,:]
			w[:,:,:]  =    W[:,:,:]
			s1[:,:,:]  =  S1[:,:,:]
			u_y[:,:,:]  =    U_y[:,:,:]
			v_y[:,:,:]  =    V_y[:,:,:]
			w_y[:,:,:]  =    W_y[:,:,:]
			s1_y[:,:,:]  =  S1_y[:,:,:]
			f.close()
		
		

#-----------------------------------------------------------------------------------------
#  YZ planes for east and west  after reading, storage format is u[t,z,y]
#-----------------------------------------------------------------------------------------
for loc in [0,1]:
	if(loc==0): root = 'east'
	if(loc==1): root = 'west'
	fname = coarse_slice_dir + root + '.nc'
	print ('Reading file ', fname)

	#---------------------------------------------------------------------------------------
	#  read in data from the desired netcdf file:  _c  means coarse quantity
	#---------------------------------------------------------------------------------------
	ncf = nc.Dataset(fname, 'r', format=fmt)
	y_c = ncf['y'][:]
	z_c = ncf['z'][:]	
	t    = ncf['t'][:]
	
	# read boundary vals
	u_c  = ncf['u'][:,:,:] 
	v_c  = ncf['v'][:,:,:]
	w_c  = ncf['w'][:,:,:]
	s1_c = ncf['s1'][:,:,:] 
	if( do_second_scalar ):
		dye_c =  ncf['s2'][:,:,:]
	
	# read boundary derivs	
	u_c_x  = ncf['u_x'][:,:,:] 
	v_c_x  = ncf['v_x'][:,:,:]
	w_c_x  = ncf['w_x'][:,:,:]
	s1_c_x = ncf['s1_x'][:,:,:] 
	if( do_second_scalar ):
		dye_c_x =  ncf['s2_x'][:,:,:]
		        
	Ly = y_c[-1]-y_c[0]                         # [m]
	Lz = z_c[-1]-z_c[0]                         # [m]
	nt = t.size
	ny_c = y_c.size
	nz_c = z_c.size
	#print 'shape of 2d data arrays u_c, etc: ', v_c.shape, nt, nz_c, ny_c  #[nt,nz_c,ny_c]        
	#print "closing ", fname
	ncf.close()
	print(' ')
	print(' ')
	
	#---------------------------------------------------------------------------------------
	#  shift origin of coarse grid and compute fine grid mesh spacing
	#---------------------------------------------------------------------------------------
	y_c = y_c - y_c[0]
	z_c = z_c - z_c[0]
	
	dy = Ly/(ny-1.)    # dy in fine grid
	dz = Lz/(nz-1.)    # dz in fine grid
	
	
	#----------------------------------------------------------------
	#  loop over the iproc and jproc processors that will write data 
	#----------------------------------------------------------------
	if( loc==0 ): iproc = 0
	if( loc==1 ): iproc = p1-1
						
	for jproc in range(p2):
			
		#-------------------------------------------------
		#  determine local fine grid coordinate vals
		#-------------------------------------------------
		y0 = 0.0 ; y1 = Ly     # y is global in YBLOCK format
		yy = np.linspace(y0,y1,ny,endpoint=True)
	
		z0 = dz*k0[jproc] ; z1 = z0 + dz*(nzvals[jproc]-1)
		zz = np.linspace(z0,z1,nzvals[jproc],endpoint=True)
		
		U  = np.zeros((nt,nzvals[jproc],ny))
		V  = np.zeros((nt,nzvals[jproc],ny))
		W  = np.zeros((nt,nzvals[jproc],ny))
		S1 = np.zeros((nt,nzvals[jproc],ny))
		
		U_x  = np.zeros((nt,nzvals[jproc],ny))
		V_x  = np.zeros((nt,nzvals[jproc],ny))
		W_x  = np.zeros((nt,nzvals[jproc],ny))
		S1_x = np.zeros((nt,nzvals[jproc],ny))
		
		#-------------------------------------------
		#  loop through time slices 
		#-------------------------------------------
		for islice in np.arange(0,nt,1):     # np.arange(0,nt,1):
	
			#-------------------------------------------
			#  compute interpolating functions over
			#  entire global coarse grid 
			#-------------------------------------------
			fu  = interpolate.interp2d(y_c, z_c,  u_c[islice,:,:].squeeze(), kind='cubic')
			fv  = interpolate.interp2d(y_c, z_c,  v_c[islice,:,:].squeeze(), kind='cubic')
			fw  = interpolate.interp2d(y_c, z_c,  w_c[islice,:,:].squeeze(), kind='cubic')
			fs1 = interpolate.interp2d(y_c, z_c, s1_c[islice,:,:].squeeze(), kind='cubic')
			
			fu_x  = interpolate.interp2d(y_c, z_c,  u_c_x[islice,:,:].squeeze(), kind='cubic')
			fv_x  = interpolate.interp2d(y_c, z_c,  v_c_x[islice,:,:].squeeze(), kind='cubic')
			fw_x  = interpolate.interp2d(y_c, z_c,  w_c_x[islice,:,:].squeeze(), kind='cubic')
			fs1_x = interpolate.interp2d(y_c, z_c, s1_c_x[islice,:,:].squeeze(), kind='cubic')
			
			#-------------------------------------------------
			# get interpolated field variables
			#-------------------------------------------------
			U[islice,:,:] = fu(yy,zz)   ; #print 'U shape ',U.shape  # ==> nt,nz,ny   
			V[islice,:,:] = fv(yy,zz) 
			W[islice,:,:] = fw(yy,zz)
			S1[islice,:,:] = fs1(yy,zz)
			
			U_x[islice,:,:] = fu_x(yy,zz)   ; #print 'U shape ',U.shape  # ==> nt,nz,ny   
			V_x[islice,:,:] = fv_x(yy,zz) 
			W_x[islice,:,:] = fw_x(yy,zz)
			S1_x[islice,:,:] = fs1_x(yy,zz)
			
		#-------------------------------------------------
		#  write local data to file
		#-------------------------------------------------				
		fname = BC_file_dir + root + '_' + str(iproc).zfill(3) + '-' + str(jproc).zfill(3) + '.nc'
		f = nc.Dataset(fname, 'w', format=fmt)
		print ('created new netcdf file: ',fname)
		timedimension = f.createDimension("timedimension", nt)
		jdimension = f.createDimension("jdimension", ny)
		kdimension = f.createDimension("kdimension", nzvals[jproc])
		time = f.createVariable('time', dtype, ("timedimension",))			
		y  = f.createVariable('y', dtype, ("jdimension",))
		z  = f.createVariable('z', dtype, ("kdimension",))
		u  = f.createVariable('u', dtype, ('timedimension','kdimension','jdimension',))
		v  = f.createVariable('v', dtype, ('timedimension','kdimension','jdimension',))
		w  = f.createVariable('w', dtype, ('timedimension','kdimension','jdimension',))
		s1 = f.createVariable('s1',dtype, ('timedimension','kdimension','jdimension',))
		if( do_second_scalar):
			s2 = f.createVariable('s2_y',dtype, ('timedimension','kdimension','jdimension',))
		u_x  = f.createVariable('u_x', dtype, ('timedimension','kdimension','jdimension',))
		v_x  = f.createVariable('v_x', dtype, ('timedimension','kdimension','jdimension',))
		w_x  = f.createVariable('w_x', dtype, ('timedimension','kdimension','jdimension',))
		s1_x = f.createVariable('s1_x',dtype, ('timedimension','kdimension','jdimension',))
		if( do_second_scalar):
			s2_x = f.createVariable('s2_x',dtype, ('timedimension','kdimension','jdimension',))
				
		time.units = 's'
		y.units = 'm'
		z.units = 'm'
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		if( do_second_scalar):
			s2.units = '1'
		u_x.units = '1/s'
		v_x.units = '1/s'
		w_x.units = '1/s'
		s1_x.units = 'kg/m4'
		if( do_second_scalar):
			s2_x.units = '1/m'
					
		# write all the data to this file
		time[:] = t[:]
		z[:] =   zz[:]
		y[:] =   yy[:]
		u[:,:,:]  =    U[:,:,:]
		v[:,:,:]  =    V[:,:,:]
		w[:,:,:]  =    W[:,:,:]
		s1[:,:,:]  =  S1[:,:,:]
		u_x[:,:,:]  =    U_x[:,:,:]
		v_x[:,:,:]  =    V_x[:,:,:]
		w_x[:,:,:]  =    W_x[:,:,:]
		s1_x[:,:,:]  =  S1_x[:,:,:]
		f.close()
			


    
print (' fine grid BC file directory              ',BC_file_dir)
print (' nested fine grid nx ny nz                ',nx,ny,nz )
print (' fine grid decomposition p1 p2            ',p1,p2,' ==> numprocs = ', p1*p2)
print (' time increment btwn BC values     [s]    ',t[1]-t[0])
print (' fine grid mesh spacing dx dy dz   [m]    ',dx,dy,dz)
print (' fine grid domain size Lx Ly Lz    [m]    ',Lx,Ly,Lz)
print (' starting/ending times in BC files [s]    ',t[0],t[-1])   


#---------------------------------------------------------
#   sanity checks afterward...
#---------------------------------------------------------

#   ls BVALS/south* | wc -l        should be p1 x p2
#   ls BVALS/north* | wc -l        should be p1 x p2  

#   ls BVALS/bottom* | wc -l        should be p1
#   ls BVALS/top* | wc -l        should be p1  

#   ls BVALS/east* | wc -l        should be p2
#   ls BVALS/west* | wc -l        should be p2
print (' '); print ('...counting the output files... '); print(' ');
cmd = 'ls  ' + BC_file_dir + 'south_* | wc -l' ; os.system(cmd)
print(' number of south_ files should be p1 x p2',p1*p2) ; print(' ')
cmd = 'ls  ' + BC_file_dir + 'north_* | wc -l' ; os.system(cmd)
print(' number of north_ files should be p1 x p2 ',p1*p2) ; print(' ')

cmd = 'ls  ' + BC_file_dir + 'bottom_* | wc -l' ; os.system(cmd)
print(' number of bottom_ files should be p1 ',p1) ; print(' ')
cmd = 'ls  ' + BC_file_dir + 'top_* | wc -l' ; os.system(cmd)
print(' number of top_ files should be p1 ',p1) ; print(' ')

cmd = 'ls  ' + BC_file_dir + 'east_* | wc -l' ; os.system(cmd)
print(' number of east_ files should be p2 ',p2) ; print(' ')
cmd = 'ls  ' + BC_file_dir + 'west_* | wc -l' ; os.system(cmd)
print(' number of west_ files should be p2 ',p2) ; print(' ')

exit()



