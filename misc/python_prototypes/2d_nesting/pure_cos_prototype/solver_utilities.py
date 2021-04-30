def create_2d_variables(nx,nz):
	import numpy as np
	a1 = np.zeros((nx,nz),dtype=float,order='F')   # leftmost index varies fastest in memory
	a2 = np.zeros((nx,nz),dtype=float,order='F')
	a3 = np.zeros((nx,nz),dtype=float,order='F')
	a4 = np.zeros((nx,nz),dtype=float,order='F')
	a5 = np.zeros((nx,nz),dtype=float,order='F')
	a6 = np.zeros((nx,nz),dtype=float,order='F')
	a7 = np.zeros((nx,nz),dtype=float,order='F')
	a8 = np.zeros((nx,nz),dtype=float,order='F')
	vars = [a1,a2,a3,a4,a5,a6,a7,a8]
	return vars

def create_1d_arrays(nx):
	import numpy as np
	a1 = np.zeros((nx,1),dtype=float)
	a2 = np.zeros((nx,1),dtype=float)
	a3 = np.zeros((nx,1),dtype=float)
	a4 = np.zeros((nx,1),dtype=float)
	a5 = np.zeros((nx,1),dtype=float)
	a6 = np.zeros((nx,1),dtype=float)
	a7 = np.zeros((nx,1),dtype=float)
	a8 = np.zeros((nx,1),dtype=float)
	vars = [a1,a2,a3,a4,a5,a6,a7,a8]
	return vars


def step(x,gamma):
#-----------------------------------------------------------------------------------------
#  Define finite width step functions near endpoints of x array
#  the step functions have zero derivative at the ends and decay from 1 at the ends
#  to 0 in the interior of the domain over the scale gamma 
#-----------------------------------------------------------------------------------------
	import numpy as np
	p = 4
	a = x[0] ; b = x[-1]
	left  = np.exp(-((x-a)/gamma)**p) 
	right = np.exp(-((x-b)/gamma)**p)
	return left,right

def ddx(f,x,frac,nvals):
#-------------------------------------------------------------------------------------------------
#    compute df/dx with f given at equally spaced grid points in [0,L] assuming good cos expansion
#
#    frac is the fraction of the high wavenumber range to filter in cos transform differentiation
#         frac<0 ==> 2/3 rule, frac0 ==> no filtering
#-------------------------------------------------------------------------------------------------
	import numpy as np
	from Fourier_stuff import dst_filtered
	nx = x.size	
	order = 1          # 1st derivative
	L = x[-1]-x[0]
	flag = -1          # even extend data, use cos series
		
	dfdx = dst_filtered(f,order,L,flag,frac)
	
	if(nvals>0):
		#-----------------------------------------------------------
		# get FD estimates at nvals pts near/including boundaries
		#-----------------------------------------------------------
		df_left,df_right = fd_deriv(f,x,nvals)		 
		#-----------------------------------------------------------
		# replace end values with finite difference estimates
		#-----------------------------------------------------------
		dfdx[0:nvals] = df_left[0:nvals]
		dfdx[nx-nvals:nx] = df_right[0:nvals]
	
		
	return dfdx


def grad(f,x,z,frac,nvals):
#-------------------------------------------------------------------------------------------------
#    compute grad f[x,z] = f_x and f_z using standard Cosine scheme f[x,z] on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(f)	
	f_x = np.zeros_like(f) 
	f_z = np.zeros_like(f)
	
	# compute f_x at each z
	for k in range(nz):
		f_x[:,k] = ddx(f[:,k],x,frac,nvals)
	
	# compute f_z at each x
	for i in range(nx):
		f_z[i,:] = ddx(f[i,:],z,frac,nvals)					
	return f_x,f_z

	
def divergence(u,w,x,z,frac,nvals):
#-------------------------------------------------------------------------------------------------
#    compute div = u_x + w_z using standard Cosine scheme 
#    u[x,z] and w[x,z]  on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(u)
	u_x = np.zeros_like(u) ; w_z = np.zeros_like(w)
	
	# compute u_x
	for k in range(nz):
		u_x[:,k] = ddx(u[:,k],x,frac,nvals)
	
	# compute w_z
	for i in range(nx):
		w_z[i,:] = ddx(w[i,:],z,frac,nvals)
	
	div = u_x + w_z	
	return div


def create_inversion_matrix(nx,nz,Lx,Lz,frac):
#-----------------------------------------------------------------------------------------
#  For solving p_xx + p_zz = f(x,z)  
#      +  homogeneous Neumann conditions on the boundary of [0,Lx] x [0,Lz] w/ (nx,nz) gridpoints
#  p and f are cosine expandable functions in both x and z directions 
#
#    p_hat = -1/(k^2+m^2) * f^hat    (except at k=m=0 where p_hat=0 )
#       ==> inversion_matrix[k,m] = -1/(k^2+m^2)
#    where this inversion is done in the Fourier space associated with the even-extended data
#
#   layout of arrays is [0 positive wavenumbers negative wavenumbers ] in each dimension
#
#-----------------------------------------------------------------------------------------
	import numpy as np
	from Fourier_stuff import fourier_wavenumbers, fourier_filter
        
	MX = (nx-1)*2 ; MZ = (nz-1)*2    # size of even extended data array, Fourier wavnumber arrays
    
	dk = np.pi/Lx
	k = fourier_wavenumbers(dk,MX)   # Fourier x wavenumbers for extended array of size MX
	kfilt = fourier_filter(k,frac)
	
       
	dm = np.pi/Lz
	m = fourier_wavenumbers(dm,MZ)   # Fourier z wavenumbers for extended array of size MZ
	mfilt = fourier_filter(m,frac)

    
	inversion_matrix = np.zeros((MX,MZ),dtype=float)
	for kk in range(MZ):
		for ii in range(MX):
			denom = (k[ii]**2 + m[kk]**2)
			if( denom==0. ):    
				inversion_matrix[ii,kk] = 0.0
			else:
				inversion_matrix[ii,kk] = (-1./denom) * kfilt[ii]*mfilt[kk]    			
	return  inversion_matrix         # return the 2d inversion matrix for later, repeated use


def poisson_invert_2d(f,inversion_matrix):
#-----------------------------------------------------------------------------------------
#  solve p_xx + p_zz = f(x,z)  
#      +  homogeneous Neumann conditions on the boundary of [0,Lx] x [0,Lz]
#  p and f are cosine expandable functions in both x and z directions 
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft2, ifft2
    from Fourier_stuff import even_extend
        
    NX,NZ = np.shape(f)                                      # size of input array
    MX = (NX-1)*2 ; MZ = (NZ-1)*2                            # size of even extended array        
    F = np.zeros((MX,MZ))                                    # array for holding 2d even extended data
    
    #-------------------------------------------------------
    # even extend data in x
    #-------------------------------------------------------
    for kk in range(NZ):
    	tmp = f[:,kk]
    	F[:,kk] = even_extend(tmp) 
    
    #-------------------------------------------------------
    # even extend data in z
    #-------------------------------------------------------
    for ii in range(MX):
    	tmp = F[ii,0:NZ]
    	F[ii,:] = even_extend(tmp)
       
    #-------------------------------------------------------
    # take FFT of 2d periodic function F
    #-------------------------------------------------------
    FHAT = fft2(F)
    
    #-----------------------------------------------------------
    # multiply elementwise by the inversion matrix -1/(k^2+m^2)
    #-----------------------------------------------------------
    p_hat = inversion_matrix * FHAT
    
    #-------------------------------------------------------------
    # take the inverse FFT to get extended, periodic version of p
    #-------------------------------------------------------------
    p = ifft2(p_hat)
    
    #-------------------------------------------------------------
    # keep and return the real portion of p in [0,Lx],[0,Lz]
    #-------------------------------------------------------------
    p = p[0:NX,0:NZ].real
    del F,FHAT,p_hat,tmp
    return p

def create_diffusion_matrix(nx,nz,Lx,Lz,frac,kappa,p,dt):
#-----------------------------------------------------------------------------------------
#  create the (extended size) diffusion array in wavenumber space given
#  the x and z diffusion coefficients kappa[0:1] and the x and z 1/2 powers p[0:1]
#
#-----------------------------------------------------------------------------------------
	import numpy as np
	from Fourier_stuff import fourier_wavenumbers, fourier_filter
        
	MX = (nx-1)*2 ; MZ = (nz-1)*2    # size of even extended data array, Fourier wavnumber arrays
    
	dk = np.pi/Lx
	k = fourier_wavenumbers(dk,MX)   # Fourier x wavenumbers for extended array of size MX
	kfilt = fourier_filter(k,frac)
	
       
	dm = np.pi/Lz
	m = fourier_wavenumbers(dm,MZ)   # Fourier z wavenumbers for extended array of size MZ
	mfilt = fourier_filter(m,frac)
    
	diffusion_matrix = np.zeros((MX,MZ),dtype=float)
	for kk in range(MZ):
		for ii in range(MX):
				xx = kappa[0]*k[ii]**(2*p[0])
				zz = kappa[1]*m[kk]**(2*p[1])
				diffusion_matrix[ii,kk] =  np.exp(-(xx+zz)*dt)   			
	return  diffusion_matrix         # return the 2d inversion matrix for later, repeated use


def diffuse_2d(f,diffusion_matrix):
#-----------------------------------------------------------------------------------------
#    exact integration for 1 time step of the diffusion term
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft2, ifft2
    from Fourier_stuff import even_extend
        
    NX,NZ = np.shape(f)                                      # size of input array
    MX = (NX-1)*2 ; MZ = (NZ-1)*2                            # size of even extended array        
    F = np.zeros((MX,MZ))                                    # array for holding 2d even extended data
    
    #-------------------------------------------------------
    # even extend data in x
    #-------------------------------------------------------
    for kk in range(NZ):
    	tmp = f[:,kk]
    	F[:,kk] = even_extend(tmp) 
    
    #-------------------------------------------------------
    # even extend data in z
    #-------------------------------------------------------
    for ii in range(MX):
    	tmp = F[ii,0:NZ]
    	F[ii,:] = even_extend(tmp)
        
    #-------------------------------------------------------
    # take FFT of 2d periodic function F
    #-------------------------------------------------------
    FHAT = fft2(F)
    
    #-----------------------------------------------------------
    # multiply elementwise by the diffusion matrix 
    #-----------------------------------------------------------
    FHAT = diffusion_matrix * FHAT
    
    #-------------------------------------------------------------
    # take the inverse FFT to get extended, periodic version of F
    #-------------------------------------------------------------
    F = ifft2(FHAT)
    
    #-------------------------------------------------------------
    # keep and return the real portion of F in [0,Lx],[0,Lz]
    #-------------------------------------------------------------
    F = F[0:NX,0:NZ].real
    return F


def boundary_smooth(f,x,z,dir,npts):
	import numpy as np
	nx = x.size  ;  nz = z.size
	dx = x[1]-x[0] ; dz = z[1]-z[0]
	
	if( dir=='x' or dir=='both' ):
		gamma = npts*dx
		for k in range(nz):
			#f[:,k] = near_boundary_smoothing(f[:,k],x,gamma)
			f[:,k] = diffuse_near_boundary(f[:,k],x,gamma)
	
	if( dir=='z' or dir=='both' ):
		gamma = npts*dz
		for i in range(nx):
			#f[i,:] = near_boundary_smoothing(f[i,:],z,gamma)
			f[i,:] = diffuse_near_boundary(f[i,:],z,gamma)
	return f
	

def near_boundary_smoothing(f,x,gamma):
#----------------------------------------------------------------------------------
#  Neighbor weights are 1 at the boundary and decrease to zero moving into
#  the interior. Well away from the boundary, the weights are zero and so there
#  is no smoothing at all. End values are left untouched.
#----------------------------------------------------------------------------------
	import numpy as np
	nx = x.size;   fs = f.copy()
	
	# define a vector of neighbor weights
	w = np.exp(-x/gamma) + np.exp(-(x[-1]-x)/gamma) 
	
	for i in range(1,nx-1):
		fs[i] = (w[i]*f[i-1] + f[i] + w[i]*f[i+1]) / (2*w[i]+1.)
	return fs	




def time_step(n,dt,U,RHS,AB_order,idx):
#  n is the current time step, time stepping starts at n=0
#  time step the eqns one step using Adams Bashforth schemes
#  toggle the indices "idx"
#  RHS[nx,nz,var_id,4]     the last index is for the 4 time levels required for AB4
	import numpy as np
	
	inc = 5
	[u,v,w,b] = U
	[MM0,MM1,MM2,MM3,method] = idx
	if( n < AB_order or np.mod(n,inc)==0 ): print("Step  ",n,method)
	
	
	if(method=='euler'):      # Euler
		u = u + dt*RHS[:,:,0,MM0]
		v = v + dt*RHS[:,:,1,MM0]
		w = w + dt*RHS[:,:,2,MM0]
		b = b + dt*RHS[:,:,3,MM0]
		if( AB_order > 1 ): method = 'AB2'
                        	
	if(n > 0 and method=='AB2'):       # AB2
		dt2 = dt/2.
		u = u + dt2*( 3.*RHS[:,:,0,MM0] - RHS[:,:,0,MM1] )
		v = v + dt2*( 3.*RHS[:,:,1,MM0] - RHS[:,:,1,MM1] )
		w = w + dt2*( 3.*RHS[:,:,2,MM0] - RHS[:,:,2,MM1] )
		b = b + dt2*( 3.*RHS[:,:,3,MM0] - RHS[:,:,3,MM1] )
		if( AB_order > 2 ): method = 'AB3'
		
	if(n > 1 and method=='AB3'):       # AB3
		dt12 = dt/12.
		u = u + dt12*( 23.*RHS[:,:,0,MM0] - 16.*RHS[:,:,0,MM1] + 5.*RHS[:,:,0,MM2] )
		v = v + dt12*( 23.*RHS[:,:,1,MM0] - 16.*RHS[:,:,1,MM1] + 5.*RHS[:,:,1,MM2] )
		w = w + dt12*( 23.*RHS[:,:,2,MM0] - 16.*RHS[:,:,2,MM1] + 5.*RHS[:,:,2,MM2] )
		b = b + dt12*( 23.*RHS[:,:,3,MM0] - 16.*RHS[:,:,3,MM1] + 5.*RHS[:,:,3,MM2] )
		if( AB_order > 3 ): method = 'AB4'
		
	if(n > 2 and method=='AB4'):        # AB4  
		dt24 = dt/24.
		u = u + dt24*( 55.*RHS[:,:,0,MM0] - 59.*RHS[:,:,0,MM1] + 37.*RHS[:,:,0,MM2] - 9.*RHS[:,:,0,MM3] )
		v = v + dt24*( 55.*RHS[:,:,1,MM0] - 59.*RHS[:,:,1,MM1] + 37.*RHS[:,:,1,MM2] - 9.*RHS[:,:,1,MM3] )
		w = w + dt24*( 55.*RHS[:,:,2,MM0] - 59.*RHS[:,:,2,MM1] + 37.*RHS[:,:,2,MM2] - 9.*RHS[:,:,2,MM3] )
		b = b + dt24*( 55.*RHS[:,:,3,MM0] - 59.*RHS[:,:,3,MM1] + 37.*RHS[:,:,3,MM2] - 9.*RHS[:,:,3,MM3] )
	
	#  Toggle M cycle pointers for AB rhs fields
	if( AB_order == 4 ):
		itmp=MM3
		MM3=MM2
		MM2=MM1
		MM1=MM0
		MM0=itmp
	elif( AB_order == 3 ):
		itmp=MM2
		MM2=MM1
		MM1=MM0
		MM0=itmp
	elif( AB_order == 2 ):
		itmp=MM1
		MM1=MM0
		MM0=itmp
	elif( AB_order == 1 ):     # flag for Euler method
		MM0 = MM0
 
	U = [u,v,w,b]                     # updated dependent variables
 	jdx = [MM0,MM1,MM2,MM3,method]    # updated index array
	return U,jdx



def zero_derivs_near_boundaries(f,x,z,npts,dir):
#-------------------------------------------------------
#  given f[x,z] w/ no particular symmetries, adjust f
#  near the boundaries so that standard cos expansions
#  can be used (and differentiated term by term) in the
#  solver algorithm
#-------------------------------------------------------
	nx = x.size ; dx=x[1]-x[0]
	nz = z.size ; dz=z[1]-z[0]
	ans = f.copy()
	
	if( dir=='x' or dir=='xz' ):
		#   adjust near x=0 and x=Lx
		gamma = npts*dx
		for k in range(nz):
			#ans[:,k] = shift_smooth_unshift(f[:,k],x,gamma)
			ans[:,k] = diffuse_near_boundary(f[:,k],x,npts)
			
	if( dir=='xz' ): f = ans        # in this case x dir already done
	
	if( dir=='z' or dir=='xz' ):
		#   adjust near z=0 and z=Lz
		gamma = npts*dz
		for i in range(nx):
			#ans[i,:] = shift_smooth_unshift(f[i,:],z,gamma)
			ans[i,:] = diffuse_near_boundary(f[i,:],z,npts)
	return ans	



def div_ustar(ustar,wstar,x,z,frac,npts,BVALS,STEP):
#-------------------------------------------------------------------------------------------------
#    compute div ustar after making necessary adjustments to accomodate straight cos scheme 
#    u[x,z] and w[x,z]  on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	
	# unpack the boundary values
	[east_vals,west_vals,bot_vals,top_vals] = BVALS
	
	[U_east,V_east,W_east,B_east] = east_vals
	[U_west,V_west,W_west,B_west] = west_vals
		
	[U_bot,V_bot,W_bot,B_bot] = bot_vals
	[U_top,V_top,W_top,B_top] = top_vals
	
	#--------------------------------------------------------------------------------------------
	# (1)  u* and w* have zero derivs at the appropriate boundaries but large divergences there
	#      that depend only on how we do the boundary adjustments, i.e. these are erroneous
	#      ===> compute the divergence, then remove the erroneous BL behavior
	#--------------------------------------------------------------------------------------------
	nvals=3*npts	# replace near bdry vals w/ FD estimates
	div = divergence(ustar,wstar,x,z,frac,nvals)
	#dir = 'xz'
	#div = extrap_near_boundaries(div,x,z,npts,dir)
	
	#-----------------------------------------------------------------------------------
	# (2) now compute the divergence associated only with the boundary "jumps"
	#     and add the two results together
	#-----------------------------------------------------------------------------------
	psi_u = construct_psi_u(x,z,[U_east,U_west,ustar],STEP)
	psi_w = construct_psi_w(x,z,[W_bot,W_top,wstar],STEP)
	
	nvals = 0	
	div = div + divergence(psi_u,psi_w,x,z,frac,nvals)
	
	#-----------------------------------------------------------------------------------
	# (3) now update the ustar vector
	#-----------------------------------------------------------------------------------
	ustar = ustar + psi_u
	wstar = wstar + psi_w		
	return div,ustar,wstar


def add_nudging(RHS,U,x,z,BVALS,STEP,dt,MM0):
	import numpy as np
	nz=z.size ; nx=x.size
	
	# unpack step functions
	[step_e,step_w,step_b,step_t] = STEP
	
	# unpack current solutions
	[u,v,w,b] = U
	
	# unpack the boundary values
	[east_vals,west_vals,bot_vals,top_vals] = BVALS
	
	[U_east,V_east,W_east,B_east] = east_vals
	[U_west,V_west,W_west,B_west] = west_vals
		
	[U_bot,V_bot,W_bot,B_bot] = bot_vals
	[U_top,V_top,W_top,B_top] = top_vals
		
	T_inv = 1./(4*dt)    # inverse of ~ fastest nudging scale
	
	for k in range(nz):
		for i in range(nx):
		
			RHS[i,k,0,MM0] = RHS[i,k,0,MM0]  \
			  - T_inv*step_e[i]*( u[i,k] - U_east[k] ) - T_inv*step_w[i]*( u[i,k] - U_west[k] ) \
			  - T_inv*step_b[k]*( u[i,k] -  U_bot[i] ) - T_inv*step_t[k]*( u[i,k] -  U_top[i] )
			  
			RHS[i,k,1,MM0] = RHS[i,k,1,MM0]  \
			  - T_inv*step_e[i]*( v[i,k] - V_east[k] ) - T_inv*step_w[i]*( v[i,k] - V_west[k] ) \
			  - T_inv*step_b[k]*( v[i,k] -  V_bot[i] ) - T_inv*step_t[k]*( v[i,k] -  V_top[i] )
			  
			RHS[i,k,2,MM0] = RHS[i,k,2,MM0]  \
			  - T_inv*step_e[i]*( w[i,k] - W_east[k] ) - T_inv*step_w[i]*( w[i,k] - W_west[k] ) \
			  - T_inv*step_b[k]*( w[i,k] -  W_bot[i] ) - T_inv*step_t[k]*( w[i,k] -  W_top[i] )
			  
			RHS[i,k,3,MM0] = RHS[i,k,3,MM0]  \
			  - T_inv*step_e[i]*( b[i,k] - B_east[k] ) - T_inv*step_w[i]*( b[i,k] - B_west[k] ) \
			  - T_inv*step_b[k]*( b[i,k] -  B_bot[i] ) - T_inv*step_t[k]*( b[i,k] -  B_top[i] )
			  			
	return RHS



def construct_psi_u(x,z,BVALS,STEP):
	import numpy as np		
	nx = x.size ; nz=z.size
	psi_u = np.zeros((nx,nz),float)
	# unpack the boundary values 
	[U_east,U_west,ustar] = BVALS
	# unpack step functions
	[step_e,step_w,step_b,step_t] = STEP
		
	#------------------------------------------------------------------------------	
	# construct an array that matches the required Dirichlet bcs at east/west
	#------------------------------------------------------------------------------
	for k in range(nz):
		a = (U_east[k]-ustar[ 0,k])
		b = (U_west[k]-ustar[-1,k])	
		psi_u[:,k] = a*step_e[:] + b*step_w[:]
	return psi_u


def construct_psi_w(x,z,BVALS,STEP):
	import numpy as np
	nx = x.size ; nz=z.size
	psi_w = np.zeros((nx,nz),float)
	# unpack the boundary values into 1d arrays and get ustar and wstar
	[W_bot,W_top,wstar] = BVALS
	# unpack step functions
	[step_e,step_w,step_b,step_t] = STEP
		
	#------------------------------------------------------------------------------	
	# construct an array that matches the required Dirichlet bcs at top/bottom
	#------------------------------------------------------------------------------
	for i in range(nx):
		a = (W_bot[i]-wstar[i, 0])
		b = (W_top[i]-wstar[i,-1])	
		psi_w[i,:] = a*step_b[:] + b*step_t[:]
	return psi_w

	
def extrap_near_boundaries(f,x,z,npts,dir):
#------------------------------------------------------------------
#  given f[x,z] with erroneous BL behavior, extrapolate
#  from interior back to boundaries to mitigate the problem, 
#  npts is the number of pts in overall near boundary adjustments
#------------------------------------------------------------------
	import numpy as np
	nx = x.size ; dx=x[1]-x[0]
	nz = z.size ; dz=z[1]-z[0]
	ans = f.copy()
	n = np.int(3*npts + 2)    # np.int(3*npts + 2) works, can it be made smaller?
	
	if( dir=='x' or dir=='xz' ):	
		for k in range(nz):
			#   adjust near x=0
			m = (f[n+1,k]-f[n,k])/dx    # slope
			x0 = x[n]
			y0 = f[n,k]
			ans[0:n,k] = y0 + m*( x[0:n]-x0 )
		
			f = ans
			#   adjust near x=Lx
			ii = nx - n - 1
			m = (f[ii+1,k]-f[ii,k])/dx    # slope
			x0 = x[ii]
			y0 = f[ii,k]
			ans[ii+1:nx,k] = y0 + m*( x[ii+1:nx]-x0 )
	
	if( dir=='xz' ): f = ans     # x direction already done...
	
	if( dir=='z' or dir=='xz' ):
		for i in range(nx):
			#   adjust near z=0
			m = (f[i,n+1]-f[i,n])/dz    # slope
			x0 = z[n]
			y0 = f[i,n]
			ans[i,0:n] = y0 + m*( z[0:n]-x0 )
		
			f = ans
			#   adjust near z=Lz
			ii = nz - n - 1
			m = (f[i,ii+1]-f[i,ii])/dz    # slope
			x0 = z[ii]
			y0 = f[i,ii]
			ans[i,ii+1:nz] = y0 + m*( z[ii+1:nz]-x0 )				
	return ans

	
def shift_smooth_unshift(f,x,gamma):
	#---------------------------------------------------------------------
	#  Slightly adjust an array of values with finite derivative at the
	#  ends to one which has zero derivs at the end. This version does
	#  NOT preserve the end values. gamma is the length scale over which
	#  the data is adjusted.
	#---------------------------------------------------------------------
	import numpy as np
	nx = x.size ; dx=x[1]-x[0]
	p = 2
	
	#------------
	# left side
	#------------
	ii = np.int((2./3.)*gamma/dx)    # trial and error value
	bval = f[ii]
	fs = f[:] - bval     # shift function to zero at left edge
	for i in range(nx):
		w = 1. - np.exp(-((x[i]-x[0])/gamma)**p)  # w 0 at edge -> 1 away
		fs[i]=fs[i]*w 
	fs[:] = fs[:] + bval
	
	f = fs
	
	#------------
	# right side 
	#------------
	ii = (nx) - np.int((2./3.)*gamma/dx)
	bval = f[ii]
	fs = f[:] - bval
	for i in range(nx):
		w = 1. - np.exp(-((x[i]-x[-1])/gamma)**p)
		fs[i]=fs[i]*w
	fs[:] = fs[:] + bval			
	return fs	


def fill_boundary_vals(x,z,t,x0,z0,WAVE_PARAMS,BVALS,BDERIVS):
	#------------------------------------------------------------
	#  Routine to store boundary values and normal derivatives
	#  in arrays that are used by the solver. Normally these
	#  values would be derived from saved parent run data but
	#  here we just evaluate the known outer solution at the
	#  boundary points.
	#------------------------------------------------------------
	import numpy as np
	nz=z.size  ; nx=x.size
	Lx = x[-1]-x[0]  ; Lz = z[-1]-z[0]
	
	# unpack so we have arrays available to fill
	[east_vals,west_vals,bot_vals,top_vals]=BVALS
	[U_east,V_east,W_east,B_east] = east_vals
	[U_west,V_west,W_west,B_west] = west_vals
	[U_bot,V_bot,W_bot,B_bot] = bot_vals
	[U_top,V_top,W_top,B_top] = top_vals
	
	[east_derivs,west_derivs,bot_derivs,top_derivs] = BDERIVS
	[U_x_east,V_x_east,W_x_east,B_x_east] = east_derivs
	[U_x_west,V_x_west,W_x_west,B_x_west] = west_derivs
	[U_z_bot,V_z_bot,W_z_bot,B_z_bot] = bot_derivs
	[U_z_top,V_z_top,W_z_top,B_z_top] = top_derivs
		
	for k in range(nz):	     #    save vals at EAST & WEST boundaries
		Z = z0+z[k]		
		id = 0           
		U_east[k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		U_west[k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		U_x_east[k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	# not actually used
		U_x_west[k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)	# not actually used
		
		id = 1           
		V_east[k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		V_west[k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		V_x_east[k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	
		V_x_west[k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)
		
		id = 2           
		W_east[k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		W_west[k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		W_x_east[k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	
		W_x_west[k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)
		
		id = 3           
		B_east[k] = parent_soln(x0   ,Z,t,id,WAVE_PARAMS)
		B_west[k] = parent_soln(x0+Lx,Z,t,id,WAVE_PARAMS)
		B_x_east[k] = parent_derivs(x0   ,Z,t,id,'x',WAVE_PARAMS)	
		B_x_west[k] = parent_derivs(x0+Lx,Z,t,id,'x',WAVE_PARAMS)	
				
		
	for i in range(nx):	    #    save vals at BOTTOM & TOP boundaries
		X = x0 + x[i]		
		# get the velocity components
		id = 0
		U_bot[i] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		U_top[i] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)
		U_z_bot[i] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		U_z_top[i] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS)
		
		id = 1
		V_bot[i] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		V_top[i] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)
		V_z_bot[i] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		V_z_top[i] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS) 
		
		id = 2
		W_bot[i] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		W_top[i] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)
		W_z_bot[i] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		W_z_top[i] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS)
		
		# get the buoyancy
		id = 3           
		B_bot[i] = parent_soln(X,z0   ,t,id,WAVE_PARAMS)
		B_top[i] = parent_soln(X,z0+Lz,t,id,WAVE_PARAMS)			
		B_z_bot[i] = parent_derivs(X,z0   ,t,id,'z',WAVE_PARAMS)
		B_z_top[i] = parent_derivs(X,z0+Lz,t,id,'z',WAVE_PARAMS)	

	# pack up the results
	east_vals=[U_east,V_east,W_east,B_east]
	west_vals=[U_west,V_west,W_west,B_west]
	bot_vals=[U_bot,V_bot,W_bot,B_bot]
	top_vals=[U_top,V_top,W_top,B_top]	
	BVALS = [east_vals,west_vals,bot_vals,top_vals]
	
	east_derivs = [U_x_east,V_x_east,W_x_east,B_x_east] 
	west_derivs = [U_x_west,V_x_west,W_x_west,B_x_west]
	bot_derivs = [U_z_bot,V_z_bot,W_z_bot,B_z_bot]
	top_derivs = [U_z_top,V_z_top,W_z_top,B_z_top] 
	BDERIVS = [east_derivs,west_derivs,bot_derivs,top_derivs]	
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

	
	
def apply_bcs(u,v,w,b,x,z,BVALS,BDERIVS,npts):
	#---------------------------------------------------------
	# In this routine, after our best attempt to get the soln,
	# we improve the behavior in the near boundary regions by
	# extending known values and derivatives into the interior
	# This is not perfect and we allow the extrapolation to
	# be shifted somewhat to match the computed solution 
	# at a point interior to the boundary layer. Boundary
	# derivatives match but the values are allowed to be slightly 
	# different than prescribed to ensure smoothness.
	# This appears to be much more robust, stable and accurate
	# than nudging towards boundary values.
	#---------------------------------------------------------	
	import numpy as np
	nx = x.size  ; nz = z.size
	
	# unpack boundary information into 1d arrays
	[east_vals,west_vals,bot_vals,top_vals]=BVALS
	[U_east,V_east,W_east,B_east] = east_vals
	[U_west,V_west,W_west,B_west] = west_vals
	[U_bot,V_bot,W_bot,B_bot] = bot_vals
	[U_top,V_top,W_top,B_top] = top_vals
	
	[east_derivs,west_derivs,bot_derivs,top_derivs] = BDERIVS
	[U_x_east,V_x_east,W_x_east,B_x_east] = east_derivs
	[U_x_west,V_x_west,W_x_west,B_x_west] = west_derivs
	[U_z_bot,V_z_bot,W_z_bot,B_z_bot] = bot_derivs
	[U_z_top,V_z_top,W_z_top,B_z_top] = top_derivs

	#------------------------------------------------------------------------
	# correct values in BL near E/W boundaries; ; match functions at x[ii]
	#------------------------------------------------------------------------
	ii=3*npts  #  6 works  i.e. 3*npts w/ npts=2
	for k in range(nz):

		# adjusting u improves u slightly but gives rise to problems in w
		f = U_east[k] + U_x_east[k]*x[:]         # low order Taylor series expansion
		offset = u[ii-1,k]-f[ii-1]
		u[0:ii,k] = f[0:ii] + offset
		
		f = U_west[k] + U_x_west[k]*(x[:]-x[-1]) # low order Taylor series expansion
		offset = u[nx-ii-1,k]-f[nx-ii-1]
		u[-ii:nx,k] = f[-ii:nx] + offset
		
		f = V_east[k] + V_x_east[k]*x[:]         # low order Taylor series expansion
		offset = v[ii-1,k]-f[ii-1]
		v[0:ii,k] = f[0:ii] + offset
		
		f = V_west[k] + V_x_west[k]*(x[:]-x[-1]) # low order Taylor series expansion
		offset = v[nx-ii-1,k]-f[nx-ii-1]
		v[-ii:nx,k] = f[-ii:nx] + offset
		
		f = W_east[k] + W_x_east[k]*x[:]         # low order Taylor series expansion
		offset = w[ii-1,k]-f[ii-1]
		w[0:ii,k] = f[0:ii] + offset
		
		f = W_west[k] + W_x_west[k]*(x[:]-x[-1]) # low order Taylor series expansion
		offset = w[nx-ii-1,k]-f[nx-ii-1]
		w[-ii:nx,k] = f[-ii:nx] + offset
		
		f = B_east[k] + B_x_east[k]*x[:]         # low order Taylor series expansion
		offset = b[ii-1,k]-f[ii-1]
		b[0:ii,k] = f[0:ii] + offset
		
		f = B_west[k] + B_x_west[k]*(x[:]-x[-1]) # low order Taylor series expansion
		offset = b[nx-ii-1,k]-f[nx-ii-1]
		b[-ii:nx,k] = f[-ii:nx] + offset
	
	#------------------------------------------------------------------------
	# correct values in BL near B/T boundaries; match functions at z[kk]
	#------------------------------------------------------------------------
	kk=3*npts #  6 works  i.e. 3*npts w/ npts=2
	for i in range(nx):
		f = U_bot[i] + U_z_bot[i]*z[:] 
		offset = u[i,kk-1] - f[kk-1]
		u[i,0:kk] = f[0:kk] + offset
		
		f = U_top[i] + U_z_top[i]*(z[:]-z[-1]) 
		offset = u[i,nz-kk-1]-f[nz-kk-1]
		u[i,-kk:nz] = f[-kk:nz] + offset
		
		f = V_bot[i] + V_z_bot[i]*z[:] 
		offset = v[i,kk-1] - f[kk-1]
		v[i,0:kk] = f[0:kk] + offset
		
		f = V_top[i] + V_z_top[i]*(z[:]-z[-1]) 
		offset = v[i,nz-kk-1]-f[nz-kk-1]
		v[i,-kk:nz] = f[-kk:nz] + offset
		
		kk=kk+2
		f = W_bot[i] + W_z_bot[i]*z[:] 
		offset = w[i,kk-1] - f[kk-1]
		w[i,0:kk] = f[0:kk] + offset
		
		f = W_top[i] + W_z_top[i]*(z[:]-z[-1]) 
		offset = w[i,nz-kk-1]-f[nz-kk-1]
		w[i,-kk:nz] = f[-kk:nz] + offset
		kk=kk-2
		
		f = B_bot[i] + B_z_bot[i]*z[:] 
		offset = b[i,kk-1] - f[kk-1]
		b[i,0:kk] = f[0:kk] + offset
		
		f = B_top[i] + B_z_top[i]*(z[:]-z[-1]) 
		offset = b[i,nz-kk-1]-f[nz-kk-1]
		b[i,-kk:nz] = f[-kk:nz] + offset
	return u,v,w,b	


def add_body_force(sigma,t,x,z,A,rhs):
	#------------------------------------------------------------
	#  add an oscillating "point source" forcing term to rhs
	#------------------------------------------------------------
	import numpy as np
	nx=x.size                       ;   nz=z.size
	dx = x[1]-x[0] ; gamma_x=2*dx   ;   dz = z[1]-z[0] ; gamma_z=2*dz  
	x0 = (x[-1]-x[0])/2.            ;   z0 = (z[-1]-z[0])/2.
	xx = A*np.cos(sigma*t)           # amplitude and time dependence
	for k in range(nz):
		for i in range(nx):
			W = np.exp(-((x[i]-x0)/gamma_x)**2) * np.exp(-((z[k]-z0)/gamma_z)**2)  # spatial shape
			rhs[i,k] = rhs[i,k] + W*xx
	return rhs
	
	
def diffuse_near_boundary(f,x,npts):
	#---------------------------------------------------------------------
	#  Slightly adjust an array of values with finite derivative at the
	#  ends to one which has zero derivs at the end. This version does
	#  NOT preserve the end values. gamma is the length scale over which
	#  the data is adjusted. Do this by mimicking diffusion. This is just
	#  a more "agnostic" version of zeroing out f_x near the boundaries.
	#  It's still O(nx) work, small compared to transforms, and so still "cheap".
	#---------------------------------------------------------------------
	import numpy as np
	nx = x.size ; h=x[1]-x[0] ; h2=h*h
	gamma = 1.5*npts*h
	nsteps = 3
	fs = f.copy()
	
	#------------------------------------------------------------------------
	# kappa=1 (but decays toward interior)
	#  integrate so that the length scale of diffusion extends to gamma
	#  do this so that explicit scheme is stable: dt = (1/3) gamma^2/kappa
	#
	#  i.e. f_t = kappa*f_xx;   integrate for
	#  a few (artificial) time steps so that smoothing extends out from
	#  the boundaries to a distance approximately equal to  gamma
	#------------------------------------------------------------------------
	kappa = np.exp(-((x-x[0])/gamma)**2) + np.exp(-((x-x[-1])/gamma)**2)
	kappa_x = -((x-x[0] )/gamma**2)*np.exp(-((x-x[0] )/gamma)**2) \
	          -((x-x[-1])/gamma**2)*np.exp(-((x-x[-1])/gamma)**2)  
			  
			  
	dt = .475*h2/np.max(kappa)   
	for n in range(nsteps):
		for i in range(1,nx-1):
			fs[i] = fs[i] + dt*kappa[i]*(fs[i-1]-2.*fs[i]+fs[i+1])/h2 +dt*kappa_x[i]*(fs[i+1]-fs[i-1])/(2.*h)
		fs[0] = fs[0] + dt*kappa[0]*(fs[1]-2.*fs[0]+fs[1])/h2                      # assume even extensions
		fs[nx-1] = fs[nx-1] + dt*kappa[nx-1]*(fs[nx-2]-2.*fs[nx-1]+fs[nx-2])/h2    # assume even extensions
	return fs	


def fd_deriv(f,x,npts):
#--------------------------------------------------------
# make FD estimates of the deriv at npts near
# and including the boundaries, save the results
# in 2  arrays of size npts
#--------------------------------------------------------
	import numpy as np
	nx = x.size
	h=x[1]-x[0]
	h12=12.*h
	df_left = np.zeros((npts,1),float)
	df_right = np.zeros((npts,1),float)
	if( npts < 3 ):
		print(" fd_deriv: npts must be >= 3 ")
		exit()
	
	# 1st 2 points
	for i in range(0,2):
		df_left[i] = (-3.*f[i] + 4.*f[i+1] -f[i+2])/(2*h)    #  error is O(h^2)
	
	# interior points on left
	for i in range(2,npts):
		df_left[i] = (f[i-2] - 8.*f[i-1] + 8.*f[i+1] - f[i+2])/h12    # error is O(h^4)
		
	# interior points on right
	ii=0
	for i in range(nx-npts,nx-2):
		df_right[ii] = (f[i-2] - 8.*f[i-1] + 8.*f[i+1] - f[i+2])/h12    # error is O(h^4)
		ii=ii+1
	
	# last 2 points	
	for i in range(nx-2,nx):
		df_right[ii] = (f[i-2] - 4.*f[i-1] + 3.*f[i])/(2*h)    #  error is O(h^2)
		ii=ii+1

	return df_left.squeeze(),df_right.squeeze()

