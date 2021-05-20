def initialize_2d_variables(nx,ny):
	import numpy as np
	a1 = np.zeros((nx,ny),float)
	a2 = np.zeros((nx,ny),float)
	a3 = np.zeros((nx,ny),float)
	a4 = np.zeros((nx,ny),float)
	a5 = np.zeros((nx,ny),float)
	a6 = np.zeros((nx,ny),float)
	a7 = np.zeros((nx,ny),float)
	a8 = np.zeros((nx,ny),float)
	return a1,a2,a3,a4,a5,a6,a7,a8
	
def fill_boundary_vals(x,y,z,t,x0,y0,z0,WAVE_PARAMS,BVALS,BDERIVS):
	#------------------------------------------------------------
	#  Routine to store boundary values and normal derivatives
	#  in arrays that are used by the solver. Normally these
	#  values would be derived from saved parent run data but
	#  here we just evaluate the known outer solution at the
	#  boundary points.
	#------------------------------------------------------------
	import numpy as np
	nz=z.size  ; nx=x.size ; ny=y.size
	Lx = x[-1]-x[0]  ; Lz = z[-1]-z[0] ; Ly=y[-1]-y[0]
	
	# unpack so we have arrays available to fill
	[east_vals,west_vals,south_vals,north_vals,bot_vals,top_vals]=BVALS
	[U_east,V_east,W_east,B_east] = east_vals
	[U_west,V_west,W_west,B_west] = west_vals	
	[U_south,V_south,W_south,B_south] = south_vals
	[U_north,V_north,W_north,B_north] = north_vals	
	[U_bot,V_bot,W_bot,B_bot] = bot_vals
	[U_top,V_top,W_top,B_top] = top_vals
	
	[east_derivs,west_derivs,south_derivs,north_derivs,bot_derivs,top_derivs] = BDERIVS
	[U_x_east,V_x_east,W_x_east,B_x_east] = east_derivs
	[U_x_west,V_x_west,W_x_west,B_x_west] = west_derivs
	[U_y_south,V_y_south,W_y_south,B_y_south] = south_derivs
	[U_y_north,V_y_north,W_y_north,B_y_north] = north_derivs
	[U_z_bot,V_z_bot,W_z_bot,B_z_bot] = bot_derivs
	[U_z_top,V_z_top,W_z_top,B_z_top] = top_derivs
		
	for k in range(nz):	     #    save vals at EAST & WEST boundaries
		Z = z0+z[k]		
		id = 0           
		U_east[:,k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		U_west[:,k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		U_x_east[:,k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	# not actually used
		U_x_west[:,k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)	# not actually used
				
		id = 1           
		V_east[:,k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		V_west[:,k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		V_x_east[:,k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	
		V_x_west[:,k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)
		
		id = 2           
		W_east[:,k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		W_west[:,k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		W_x_east[:,k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	
		W_x_west[:,k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)
		
		id = 3           
		B_east[:,k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		B_west[:,k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		B_x_east[:,k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	
		B_x_west[:,k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)	
	
	for i in range(nx):	     #    save vals at SOUTH & NORTH boundaries
		X = x0+x[i]
		for k in range(nz):
			Z = z0+z[k]
							
			id = 0           
			U_south[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			U_north[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			U_y_south[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)	# not actually used
			U_y_north[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)	# not actually used
		
			id = 1           
			V_south[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			V_north[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			V_y_south[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)	
			V_y_north[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)
		
			id = 2           
			W_south[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			W_north[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			W_y_south[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)	
			W_y_north[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)
		
			id = 3           
			B_south[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			B_north[i,k] = parent_soln(X,Z,t,id,WAVE_PARAMS)
			B_y_south[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)	
			B_y_north[i,k] = parent_derivs(X,Z,t,id,'y',WAVE_PARAMS)
					
		
	for i in range(nx):	    #    save vals at BOTTOM & TOP boundaries
		X = x0 + x[i]		
		# get the velocity components
		id = 0
		U_bot[i,:] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		U_top[i,:] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)
		U_z_bot[i,:] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		U_z_top[i,:] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS)
		
		id = 1
		V_bot[i,:] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		V_top[i,:] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)
		V_z_bot[i,:] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		V_z_top[i,:] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS) 
		
		id = 2
		W_bot[i,:] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		W_top[i,:] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)
		W_z_bot[i,:] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		W_z_top[i,:] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS)
		
		# get the buoyancy
		id = 3           
		B_bot[i,:] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		B_top[i,:] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)			
		B_z_bot[i,:] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		B_z_top[i,:] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS)	

	# pack up the results
	east_vals=[U_east,V_east,W_east,B_east]
	west_vals=[U_west,V_west,W_west,B_west]
	south_vals=[U_south,V_south,W_south,B_south]
	north_vals=[U_north,V_north,W_north,B_north]
	bot_vals=[U_bot,V_bot,W_bot,B_bot]
	top_vals=[U_top,V_top,W_top,B_top]	
	BVALS = [east_vals,west_vals,south_vals,north_vals,bot_vals,top_vals]
	
	east_derivs = [U_x_east,V_x_east,W_x_east,B_x_east] 
	west_derivs = [U_x_west,V_x_west,W_x_west,B_x_west]
	south_derivs = [U_y_south,V_y_south,W_y_south,B_y_south]
	north_derivs = [U_y_north,V_y_north,W_y_north,B_y_north]
	bot_derivs = [U_z_bot,V_z_bot,W_z_bot,B_z_bot]
	top_derivs = [U_z_top,V_z_top,W_z_top,B_z_top] 
	BDERIVS = [east_derivs,west_derivs,south_derivs,north_derivs,bot_derivs,top_derivs]	
	return BVALS,BDERIVS
	
def parent_soln(X,Z,t,id,WAVE_PARAMS):
	#-----------------------------------------------------------	
	# function to evaluate parent soln at X,Z,t  ; 
	# here an IW mode with wave parameters available globally
	# id = [0,1,2,3,6] for [U,V,W,B,P] respectively
	#-----------------------------------------------------------
	import numpy as np
	[A,f,N2,omega,k_iw,m_iw,phase] = WAVE_PARAMS
	f2 = f*f ; omega2 = omega*omega
	argz = m_iw*Z
	argx = k_iw*X - omega*t + phase
	if(id==0): ans = A*np.cos(argz)*np.cos(argx) # U
	if(id==1): ans = A*(f/omega)*np.cos(argz)*np.sin(argx) # V    
	if(id==2): ans = A*(k_iw/m_iw)*np.sin(argz)*np.sin(argx) # W 
	if(id==3): ans = -A*(k_iw/m_iw)*(N2/omega)*np.sin(argz)*np.cos(argx) # B
	if(id==6): ans = -A*(1./(k_iw*omega))*(f2-omega2)*np.cos(argz)*np.cos(argx) # P	
	return ans

def parent_derivs(X,Z,t,id,dir,WAVE_PARAMS):
	#----------------------------------------------------------------------
	# function to evaluate certain derivs of parent soln at X,Z,t 
	# id = 0,1,2,3 for U,V,W,B
	# dir = 'x','y' or 'z'
	# used for normal derivs of tangential vels at boundaries
	#----------------------------------------------------------------------
	import numpy as np
	[A,f,N2,omega,k_iw,m_iw,phase] = WAVE_PARAMS
	f2 = f*f ; omega2 = omega*omega
	argz = m_iw*Z
	argx = k_iw*X - omega*t + phase	
	
	if(id==0 and dir=='x'):    
		ans = -k_iw*A*np.cos(argz)*np.sin(argx)            # U_x
	if(id==0 and dir=='z'):    
		ans = -m_iw*A*np.sin(argz)*np.cos(argx)            # U_z
		
	if(id==1 and dir=='x'):    
		ans = k_iw*A*(f/omega)*np.cos(argz)*np.cos(argx)    # V_x
	if(id==1 and dir=='z'):    
		ans = -m_iw*A*(f/omega)*np.sin(argz)*np.sin(argx)   # V_z
		
	if(id==2 and dir=='x'):    
		ans = k_iw*A*(k_iw/m_iw)*np.sin(argz)*np.cos(argx)  # W_x
	if(id==2 and dir=='z'):    
		ans = m_iw*A*(k_iw/m_iw)*np.cos(argz)*np.sin(argx)  # W_z
		
	if(id==3 and dir=='x'):    
		ans = k_iw*A*(k_iw/m_iw)*(N2/omega)*np.sin(argz)*np.sin(argx)  # B_x
	if(id==3 and dir=='z'):    
		ans = -m_iw*A*(k_iw/m_iw)*(N2/omega)*np.cos(argz)*np.cos(argx)  # B_z
		
	if(dir=='y'):    ans = 0. 
	return ans

def write_east_vals(tval,islice,yvec,zvec,east_vals,east_derivs,outdir,nt):
	#---------------------------------------------------------------------------------------
	#  import and rename the various modules to be used
	#---------------------------------------------------------------------------------------
	import os,sys
	import numpy as np
	import netCDF4 as nc
	
	fmt = "NETCDF4"
	dtype = np.float64
	do_second_scalar=False
	 	
	ny = yvec.size
	nz = zvec.size
	
	fname = outdir + 'east.nc' 
	if(islice==0): 
		mode = 'w'             
	else:
		mode = 'a'
    
    # open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(fname, mode, format=fmt) 
        
	if( islice==0 ):        
		tdim = nc.createDimension("tdim", nt)    # None or 0 for unlimited dimension
		jdim = nc.createDimension("jdim", ny)
		kdim = nc.createDimension("kdim", nz)
	
		# A variable has a name, a type, a shape and some data values. The shape of the variable can be stated 	
		# using the tuple of the dimension names. Can also add attributes.
		t = nc.createVariable('t', dtype, ("tdim",))
		y = nc.createVariable('y', dtype, ("jdim",))
		z = nc.createVariable('z', dtype, ("kdim",))
		t.units = "s"
		y.units = "m"
		z.units = "m"
		
		# Write the dimension data
		y[:] = yvec[:]
		z[:] = zvec[:]
        
		u  = nc.createVariable('u', dtype, ('tdim','kdim','jdim',))
		v  = nc.createVariable('v', dtype, ('tdim','kdim','jdim',))
		w  = nc.createVariable('w', dtype, ('tdim','kdim','jdim',))
		s1 = nc.createVariable('s1',dtype, ('tdim','kdim','jdim',))
		if( do_second_scalar):
			s2 = f.createVariable('s2','d', ('tdim','kdim','jdim',))
		#print ("-- Created variables")
		
		u_x  = nc.createVariable('u_x', dtype, ('tdim','kdim','jdim',))
		v_x  = nc.createVariable('v_x', dtype, ('tdim','kdim','jdim',))
		w_x  = nc.createVariable('w_x', dtype, ('tdim','kdim','jdim',))
		s1_x = nc.createVariable('s1_x',dtype, ('tdim','kdim','jdim',))
		if( do_second_scalar):
			s2_x = f.createVariable('s2_x','d', ('tdim','kdim','jdim',))
		#print ("-- Created variables")
    
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		u_x.units = '1/s'
		v_x.units = '1/s'
		w_x.units = '1/s'
		s1_x.units = 'kg/m4'
		if( do_second_scalar):
			s1.units = '1'
			s1_x.units = '1/m'
		#print ("-- attached units strings") 
	else:
		t   = nc.variables['t']
		u   = nc.variables['u']
		v   = nc.variables['v']
		w   = nc.variables['w']
		s1  = nc.variables['s1']
		u_x = nc.variables['u_x']
		v_x = nc.variables['v_x']
		w_x = nc.variables['w_x']
		s1_x = nc.variables['s1_x']
    
    # unpack the east boundary values    
	[U_east,V_east,W_east,B_east] = east_vals
	
	# unpack the east normal derivatives
	[U_x_east,V_x_east,W_x_east,B_x_east] = east_derivs
	
	g=9.81; rho0=1027.  
	S1_east = -(rho0/g) * B_east 
	S1_x_east = -(rho0/g) * B_x_east   
               
	#----------------------------------------------------------------------
	#  write data into fname with the time slice equal to islice
	#---------------------------------------------------------------------- 
	t[islice] = tval   
	u[islice,:,:]    =  U_east.transpose()          # U_east(y,z) ==> (z,y)
	v[islice,:,:]    =  V_east.transpose()
	w[islice,:,:]    =  W_east.transpose()
	s1[islice,:,:]   =  S1_east.transpose()
	u_x[islice,:,:]  =  U_x_east.transpose()      # U_east(y,z) ==> (z,y)
	v_x[islice,:,:]  =  V_x_east.transpose()
	w_x[islice,:,:]  =  W_x_east.transpose()
	s1_x[islice,:,:] =  S1_x_east.transpose()
	if( do_second_scalar):
		s2[islice,:,:,:]   = S2_east.transpose()
		s2_x[islice,:,:,:] = S2_x_east.transpose()
	
    
	#print("-- filled/wrote east_vals and east_derivs for islice ",islice)
	nc.close()
	success = True	
	return success

def write_west_vals(tval,islice,yvec,zvec,west_vals,west_derivs,outdir,nt):
	#---------------------------------------------------------------------------------------
	#  import and rename the various modules to be used
	#---------------------------------------------------------------------------------------
	import os,sys
	import numpy as np
	import netCDF4 as nc
	
	fmt = "NETCDF4"
	dtype = np.float64
	do_second_scalar=False
	 	
	ny = yvec.size
	nz = zvec.size
	
	fname = outdir + 'west.nc' 
	if(islice==0): 
		mode = 'w'             
	else:
		mode = 'a'
    
    # open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(fname, mode, format=fmt) 
        
	if( islice==0 ):        
		tdim = nc.createDimension("tdim", nt)    # None or 0 for unlimited dimension
		jdim = nc.createDimension("jdim", ny)
		kdim = nc.createDimension("kdim", nz)
	
		# A variable has a name, a type, a shape and some data values. The shape of the variable can be stated 	
		# using the tuple of the dimension names. Can also add attributes.
		t = nc.createVariable('t', dtype, ("tdim",))
		y = nc.createVariable('y', dtype, ("jdim",))
		z = nc.createVariable('z', dtype, ("kdim",))
		t.units = "s"
		y.units = "m"
		z.units = "m"
		
		# Write the dimension data
		y[:] = yvec[:]
		z[:] = zvec[:]
        
		u  = nc.createVariable('u', dtype, ('tdim','kdim','jdim',))
		v  = nc.createVariable('v', dtype, ('tdim','kdim','jdim',))
		w  = nc.createVariable('w', dtype, ('tdim','kdim','jdim',))
		s1 = nc.createVariable('s1',dtype, ('tdim','kdim','jdim',))
		if( do_second_scalar):
			s2 = f.createVariable('s2','d', ('tdim','kdim','jdim',))
		#print ("-- Created variables")
		
		u_x  = nc.createVariable('u_x', dtype, ('tdim','kdim','jdim',))
		v_x  = nc.createVariable('v_x', dtype, ('tdim','kdim','jdim',))
		w_x  = nc.createVariable('w_x', dtype, ('tdim','kdim','jdim',))
		s1_x = nc.createVariable('s1_x',dtype, ('tdim','kdim','jdim',))
		if( do_second_scalar):
			s2_x = f.createVariable('s2_x','d', ('tdim','kdim','jdim',))
		#print ("-- Created variables")
    
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		u_x.units = '1/s'
		v_x.units = '1/s'
		w_x.units = '1/s'
		s1_x.units = 'kg/m4'
		if( do_second_scalar):
			s1.units = '1'
			s1_x.units = '1/m'
		#print ("-- attached units strings") 
	else:
		t   = nc.variables['t']
		u   = nc.variables['u']
		v   = nc.variables['v']
		w   = nc.variables['w']
		s1  = nc.variables['s1']
		u_x = nc.variables['u_x']
		v_x = nc.variables['v_x']
		w_x = nc.variables['w_x']
		s1_x = nc.variables['s1_x']
    
    # unpack the west boundary values    
	[U_west,V_west,W_west,B_west] = west_vals
	
	# unpack the west normal derivatives
	[U_x_west,V_x_west,W_x_west,B_x_west] = west_derivs
	
	g=9.81; rho0=1027.  
	S1_west = -(rho0/g) * B_west 
	S1_x_west = -(rho0/g) * B_x_west   
               
	#----------------------------------------------------------------------
	#  write data into fname with the time slice equal to islice
	#---------------------------------------------------------------------- 
	t[islice] = tval   
	u[islice,:,:]    =  U_west.transpose()          # U_west(y,z) ==> (z,y)
	v[islice,:,:]    =  V_west.transpose()
	w[islice,:,:]    =  W_west.transpose()
	s1[islice,:,:]   =  S1_west.transpose()
	u_x[islice,:,:]  =  U_x_west.transpose()      # U_west(y,z) ==> (z,y)
	v_x[islice,:,:]  =  V_x_west.transpose()
	w_x[islice,:,:]  =  W_x_west.transpose()
	s1_x[islice,:,:] =  S1_x_west.transpose()
	if( do_second_scalar):
		s2[islice,:,:,:]   = S2_west.transpose()
		s2_x[islice,:,:,:] = S2_x_west.transpose()
	
    
	#print("-- filled/wrote west_vals and west_derivs for islice ",islice)
	nc.close()
	success = True	
	return success


def write_south_vals(tval,islice,xvec,zvec,south_vals,south_derivs,outdir,nt):
	#---------------------------------------------------------------------------------------
	#  import and rename the various modules to be used
	#---------------------------------------------------------------------------------------
	import os,sys
	import numpy as np
	import netCDF4 as nc
	
	fmt = "NETCDF4"
	dtype = np.float64
	do_second_scalar=False
	 	
	nx = xvec.size
	nz = zvec.size
	
	fname = outdir + 'south.nc' 
	if(islice==0): 
		mode = 'w'             
	else:
		mode = 'a'
    
    # open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(fname, mode, format=fmt) 
        
	if( islice==0 ):        
		tdim = nc.createDimension("tdim", nt)    # None or 0 for unlimited dimension
		idim = nc.createDimension("idim", nx)
		kdim = nc.createDimension("kdim", nz)
	
		# A variable has a name, a type, a shape and some data values. The shape of the variable can be stated 	
		# using the tuple of the dimension names. Can also add attributes.
		t = nc.createVariable('t', dtype, ("tdim",))
		x = nc.createVariable('x', dtype, ("idim",))
		z = nc.createVariable('z', dtype, ("kdim",))
		t.units = "s"
		x.units = "m"
		z.units = "m"
		
		# Write the dimension data
		x[:] = xvec[:]
		z[:] = zvec[:]
        
		u  = nc.createVariable('u', dtype, ('tdim','kdim','idim',))
		v  = nc.createVariable('v', dtype, ('tdim','kdim','idim',))
		w  = nc.createVariable('w', dtype, ('tdim','kdim','idim',))
		s1 = nc.createVariable('s1',dtype, ('tdim','kdim','idim',))
		if( do_second_scalar):
			s2 = f.createVariable('s2','d', ('tdim','kdim','idim',))
		#print ("-- Created variables")
		
		u_y  = nc.createVariable('u_y', dtype, ('tdim','kdim','idim',))
		v_y  = nc.createVariable('v_y', dtype, ('tdim','kdim','idim',))
		w_y  = nc.createVariable('w_y', dtype, ('tdim','kdim','idim',))
		s1_y = nc.createVariable('s1_y',dtype, ('tdim','kdim','idim',))
		if( do_second_scalar):
			s2_y = f.createVariable('s2_y','d', ('tdim','kdim','idim',))
		#print ("-- Created variables")
    
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		u_y.units = '1/s'
		v_y.units = '1/s'
		w_y.units = '1/s'
		s1_y.units = 'kg/m4'
		if( do_second_scalar):
			s1.units = '1'
			s1_y.units = '1/m'
		#print ("-- attached units strings") 
	else:
		t   = nc.variables['t']
		u   = nc.variables['u']
		v   = nc.variables['v']
		w   = nc.variables['w']
		s1  = nc.variables['s1']
		u_y = nc.variables['u_y']
		v_y = nc.variables['v_y']
		w_y = nc.variables['w_y']
		s1_y = nc.variables['s1_y']
    
    # unpack the south boundary values    
	[U_south,V_south,W_south,B_south] = south_vals
	
	# unpack the south normal derivatives
	[U_y_south,V_y_south,W_y_south,B_y_south] = south_derivs
	
	g=9.81; rho0=1027.  
	S1_south = -(rho0/g) * B_south 
	S1_y_south = -(rho0/g) * B_y_south   
               
	#----------------------------------------------------------------------
	#  write data into fname with the time slice equal to islice
	#---------------------------------------------------------------------- 
	t[islice] = tval   
	u[islice,:,:]    =  U_south.transpose()          # U_south(yx,z) ==> (z,x)
	v[islice,:,:]    =  V_south.transpose()
	w[islice,:,:]    =  W_south.transpose()
	s1[islice,:,:]   =  S1_south.transpose()
	u_y[islice,:,:]  =  U_y_south.transpose()      # U_south(x,z) ==> (z,x)
	v_y[islice,:,:]  =  V_y_south.transpose()
	w_y[islice,:,:]  =  W_y_south.transpose()
	s1_y[islice,:,:] =  S1_y_south.transpose()
	if( do_second_scalar):
		s2[islice,:,:,:]   = S2_south.transpose()
		s2_y[islice,:,:,:] = S2_y_south.transpose()
	
    
	#print("-- filled/wrote south_vals and south_derivs for islice ",islice)
	nc.close()
	success = True	
	return success
	
def write_north_vals(tval,islice,xvec,zvec,north_vals,north_derivs,outdir,nt):
	#---------------------------------------------------------------------------------------
	#  import and rename the various modules to be used
	#---------------------------------------------------------------------------------------
	import os,sys
	import numpy as np
	import netCDF4 as nc
	
	fmt = "NETCDF4"
	dtype = np.float64
	do_second_scalar=False
	 	
	nx = xvec.size
	nz = zvec.size
	
	fname = outdir + 'north.nc' 
	if(islice==0): 
		mode = 'w'             
	else:
		mode = 'a'
    
    # open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(fname, mode, format=fmt) 
        
	if( islice==0 ):        
		tdim = nc.createDimension("tdim", nt)    # None or 0 for unlimited dimension
		idim = nc.createDimension("idim", nx)
		kdim = nc.createDimension("kdim", nz)
	
		# A variable has a name, a type, a shape and some data values. The shape of the variable can be stated 	
		# using the tuple of the dimension names. Can also add attributes.
		t = nc.createVariable('t', dtype, ("tdim",))
		x = nc.createVariable('x', dtype, ("idim",))
		z = nc.createVariable('z', dtype, ("kdim",))
		t.units = "s"
		x.units = "m"
		z.units = "m"
		
		# Write the dimension data
		x[:] = xvec[:]
		z[:] = zvec[:]
        
		u  = nc.createVariable('u', dtype, ('tdim','kdim','idim',))
		v  = nc.createVariable('v', dtype, ('tdim','kdim','idim',))
		w  = nc.createVariable('w', dtype, ('tdim','kdim','idim',))
		s1 = nc.createVariable('s1',dtype, ('tdim','kdim','idim',))
		if( do_second_scalar):
			s2 = f.createVariable('s2','d', ('tdim','kdim','idim',))
		#print ("-- Created variables")
		
		u_y  = nc.createVariable('u_y', dtype, ('tdim','kdim','idim',))
		v_y  = nc.createVariable('v_y', dtype, ('tdim','kdim','idim',))
		w_y  = nc.createVariable('w_y', dtype, ('tdim','kdim','idim',))
		s1_y = nc.createVariable('s1_y',dtype, ('tdim','kdim','idim',))
		if( do_second_scalar):
			s2_y = f.createVariable('s2_y','d', ('tdim','kdim','idim',))
		#print ("-- Created variables")
    
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		u_y.units = '1/s'
		v_y.units = '1/s'
		w_y.units = '1/s'
		s1_y.units = 'kg/m4'
		if( do_second_scalar):
			s1.units = '1'
			s1_y.units = '1/m'
		#print ("-- attached units strings") 
	else:
		t   = nc.variables['t']
		u   = nc.variables['u']
		v   = nc.variables['v']
		w   = nc.variables['w']
		s1  = nc.variables['s1']
		u_y = nc.variables['u_y']
		v_y = nc.variables['v_y']
		w_y = nc.variables['w_y']
		s1_y = nc.variables['s1_y']
    
    # unpack the north boundary values    
	[U_north,V_north,W_north,B_north] = north_vals
	
	# unpack the north normal derivatives
	[U_y_north,V_y_north,W_y_north,B_y_north] = north_derivs
	
	g=9.81; rho0=1027.  
	S1_north = -(rho0/g) * B_north 
	S1_y_north = -(rho0/g) * B_y_north   
               
	#----------------------------------------------------------------------
	#  write data into fname with the time slice equal to islice
	#---------------------------------------------------------------------- 
	t[islice] = tval   
	u[islice,:,:]    =  U_north.transpose()          # U_north(yx,z) ==> (z,x)
	v[islice,:,:]    =  V_north.transpose()
	w[islice,:,:]    =  W_north.transpose()
	s1[islice,:,:]   =  S1_north.transpose()
	u_y[islice,:,:]  =  U_y_north.transpose()      # U_north(x,z) ==> (z,x)
	v_y[islice,:,:]  =  V_y_north.transpose()
	w_y[islice,:,:]  =  W_y_north.transpose()
	s1_y[islice,:,:] =  S1_y_north.transpose()
	if( do_second_scalar):
		s2[islice,:,:,:]   = S2_north.transpose()
		s2_y[islice,:,:,:] = S2_y_north.transpose()
	
    
	#print("-- filled/wrote north_vals and north_derivs for islice ",islice)
	nc.close()
	success = True	
	return success


def write_bottom_vals(tval,islice,xvec,yvec,bottom_vals,bottom_derivs,outdir,nt):
	#---------------------------------------------------------------------------------------
	#  import and rename the various modules to be used
	#---------------------------------------------------------------------------------------
	import os,sys
	import numpy as np
	import netCDF4 as nc
	
	fmt = "NETCDF4"
	dtype = np.float64
	do_second_scalar=False
	 	
	nx = xvec.size
	ny = yvec.size
	
	fname = outdir + 'bottom.nc' 
	if(islice==0): 
		mode = 'w'             
	else:
		mode = 'a'
    
    # open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(fname, mode, format=fmt) 
        
	if( islice==0 ):        
		tdim = nc.createDimension("tdim", nt)    # None or 0 for unlimited dimension
		idim = nc.createDimension("idim", nx)
		jdim = nc.createDimension("jdim", ny)
	
		# A variable has a name, a type, a shape and some data values. The shape of the variable can be stated 	
		# using the tuple of the dimension names. Can also add attributes.
		t = nc.createVariable('t', dtype, ("tdim",))
		x = nc.createVariable('x', dtype, ("idim",))
		y = nc.createVariable('y', dtype, ("jdim",))
		t.units = "s"
		x.units = "m"
		y.units = "m"
		
		# Write the dimension data
		x[:] = xvec[:]
		y[:] = yvec[:]
        
		u  = nc.createVariable('u', dtype, ('tdim','jdim','idim',))
		v  = nc.createVariable('v', dtype, ('tdim','jdim','idim',))
		w  = nc.createVariable('w', dtype, ('tdim','jdim','idim',))
		s1 = nc.createVariable('s1',dtype, ('tdim','jdim','idim',))
		if( do_second_scalar):
			s2 = f.createVariable('s2','d', ('tdim','jdim','idim',))
		#print ("-- Created variables")
		
		u_z  = nc.createVariable('u_z', dtype, ('tdim','jdim','idim',))
		v_z  = nc.createVariable('v_z', dtype, ('tdim','jdim','idim',))
		w_z  = nc.createVariable('w_z', dtype, ('tdim','jdim','idim',))
		s1_z = nc.createVariable('s1_z',dtype, ('tdim','jdim','idim',))
		if( do_second_scalar):
			s2_z = f.createVariable('s2_z','d', ('tdim','jdim','idim',))
		#print ("-- Created variables")
    
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		u_z.units = '1/s'
		v_z.units = '1/s'
		w_z.units = '1/s'
		s1_z.units = 'kg/m4'
		if( do_second_scalar):
			s1.units = '1'
			s1_z.units = '1/m'
		#print ("-- attached units strings") 
	else:
		t   = nc.variables['t']
		u   = nc.variables['u']
		v   = nc.variables['v']
		w   = nc.variables['w']
		s1  = nc.variables['s1']
		u_z = nc.variables['u_z']
		v_z = nc.variables['v_z']
		w_z = nc.variables['w_z']
		s1_z = nc.variables['s1_z']
    
    # unpack the bottom boundary values    
	[U_bottom,V_bottom,W_bottom,B_bottom] = bottom_vals
	
	# unpack the bottom normal derivatives
	[U_z_bottom,V_z_bottom,W_z_bottom,B_z_bottom] = bottom_derivs
	
	g=9.81; rho0=1027.  
	S1_bottom = -(rho0/g) * B_bottom 
	S1_z_bottom = -(rho0/g) * B_z_bottom   
               
	#----------------------------------------------------------------------
	#  write data into fname with the time slice equal to islice
	#---------------------------------------------------------------------- 
	t[islice] = tval   
	u[islice,:,:]    =  U_bottom.transpose()          # U_bottom(x,y) ==> (y,x)
	v[islice,:,:]    =  V_bottom.transpose()
	w[islice,:,:]    =  W_bottom.transpose()
	s1[islice,:,:]   =  S1_bottom.transpose()
	u_z[islice,:,:]  =  U_z_bottom.transpose()        # U_bottom(x,y) ==> (y,x)
	v_z[islice,:,:]  =  V_z_bottom.transpose()
	w_z[islice,:,:]  =  W_z_bottom.transpose()
	s1_z[islice,:,:] =  S1_z_bottom.transpose()
	if( do_second_scalar):
		s2[islice,:,:,:]   = S2_bottom.transpose()
		s2_z[islice,:,:,:] = S2_z_bottom.transpose()
	
    
	#print("-- filled/wrote bottom_vals and bottom_derivs for islice ",islice)
	nc.close()
	success = True	
	return success

def write_top_vals(tval,islice,xvec,yvec,top_vals,top_derivs,outdir,nt):
	#---------------------------------------------------------------------------------------
	#  import and rename the various modules to be used
	#---------------------------------------------------------------------------------------
	import os,sys
	import numpy as np
	import netCDF4 as nc
	
	fmt = "NETCDF4"
	dtype = np.float64
	do_second_scalar=False
	 	
	nx = xvec.size
	ny = yvec.size
	
	fname = outdir + 'top.nc' 
	if(islice==0): 
		mode = 'w'             
	else:
		mode = 'a'
    
    # open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(fname, mode, format=fmt) 
        
	if( islice==0 ):        
		tdim = nc.createDimension("tdim", nt)    # None or 0 for unlimited dimension
		idim = nc.createDimension("idim", nx)
		jdim = nc.createDimension("jdim", ny)
	
		# A variable has a name, a type, a shape and some data values. The shape of the variable can be stated 	
		# using the tuple of the dimension names. Can also add attributes.
		t = nc.createVariable('t', dtype, ("tdim",))
		x = nc.createVariable('x', dtype, ("idim",))
		y = nc.createVariable('y', dtype, ("jdim",))
		t.units = "s"
		x.units = "m"
		y.units = "m"
		
		# Write the dimension data
		x[:] = xvec[:]
		y[:] = yvec[:]
        
		u  = nc.createVariable('u', dtype, ('tdim','jdim','idim',))
		v  = nc.createVariable('v', dtype, ('tdim','jdim','idim',))
		w  = nc.createVariable('w', dtype, ('tdim','jdim','idim',))
		s1 = nc.createVariable('s1',dtype, ('tdim','jdim','idim',))
		if( do_second_scalar):
			s2 = f.createVariable('s2','d', ('tdim','jdim','idim',))
		#print ("-- Created variables")
		
		u_z  = nc.createVariable('u_z', dtype, ('tdim','jdim','idim',))
		v_z  = nc.createVariable('v_z', dtype, ('tdim','jdim','idim',))
		w_z  = nc.createVariable('w_z', dtype, ('tdim','jdim','idim',))
		s1_z = nc.createVariable('s1_z',dtype, ('tdim','jdim','idim',))
		if( do_second_scalar):
			s2_z = f.createVariable('s2_z','d', ('tdim','jdim','idim',))
		#print ("-- Created variables")
    
		u.units = 'm/s'
		v.units = 'm/s'
		w.units = 'm/s'
		s1.units = 'kg/m3'
		u_z.units = '1/s'
		v_z.units = '1/s'
		w_z.units = '1/s'
		s1_z.units = 'kg/m4'
		if( do_second_scalar):
			s1.units = '1'
			s1_z.units = '1/m'
		#print ("-- attached units strings") 
	else:
		t   = nc.variables['t']
		u   = nc.variables['u']
		v   = nc.variables['v']
		w   = nc.variables['w']
		s1  = nc.variables['s1']
		u_z = nc.variables['u_z']
		v_z = nc.variables['v_z']
		w_z = nc.variables['w_z']
		s1_z = nc.variables['s1_z']
    
    # unpack the top boundary values    
	[U_top,V_top,W_top,B_top] = top_vals
	
	# unpack the top normal derivatives
	[U_z_top,V_z_top,W_z_top,B_z_top] = top_derivs
	
	g=9.81; rho0=1027.  
	S1_top = -(rho0/g) * B_top 
	S1_z_top = -(rho0/g) * B_z_top   
               
	#----------------------------------------------------------------------
	#  write data into fname with the time slice equal to islice
	#---------------------------------------------------------------------- 
	t[islice] = tval   
	u[islice,:,:]    =  U_top.transpose()          # U_top(x,y) ==> (y,x)
	v[islice,:,:]    =  V_top.transpose()
	w[islice,:,:]    =  W_top.transpose()
	s1[islice,:,:]   =  S1_top.transpose()
	u_z[islice,:,:]  =  U_z_top.transpose()        # U_top(x,y) ==> (y,x)
	v_z[islice,:,:]  =  V_z_top.transpose()
	w_z[islice,:,:]  =  W_z_top.transpose()
	s1_z[islice,:,:] =  S1_z_top.transpose()
	if( do_second_scalar):
		s2[islice,:,:,:]   = S2_top.transpose()
		s2_z[islice,:,:,:] = S2_z_top.transpose()
	
    
	#print("-- filled/wrote top_vals and top_derivs for islice ",islice)
	nc.close()
	success = True	
	return success
