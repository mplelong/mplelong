def prepare_BernoulliCosine(x,Q,inc):
#-------------------------------------------------------------------------------------------------
#    construct and factor the matrices required for the Bernoulli expansion coefficients
#    LU = [lu_A,lu_B,piv_A,piv_B] the prefactored matrices constraining the expansion coeffs
#-------------------------------------------------------------------------------------------------
	import numpy as np
	from Bernoulli_polynomials import setup_factor_matrix
	x0 = x[0]     # even extension of f has discontinuous deriv at x=0
	lu_A,piv_A = setup_factor_matrix(x,Q,x0,inc)
	x0 = x[-1]    # even extension of f also has discontinuous deriv at x=L
	lu_B,piv_B = setup_factor_matrix(x,Q,x0,inc)
	LU = [lu_A,lu_B,piv_A,piv_B]
	return LU

def evaluate_basis_functions(x,Q,k):
#-------------------------------------------------------------------------------------------------
#    Construct the U_n basis functions for n=1,3,5... and their derivatives
#    at each x point for the (Q+1)/2 term series
#    The S_0 series are expended about x=0 while the S_L series are expanded about x=x[-1]=L
#-------------------------------------------------------------------------------------------------
	import numpy as np
	from Bernoulli_polynomials import U_n,ddx_U_n,kth_deriv_U_n
	nx = x.size
	# number of terms in the series expansions
	M = (Q+1)/2
	U_0 = np.zeros( (nx,M),float )   ; U_L = np.zeros( (nx,M),float )
	U_0_x = np.zeros( (nx,M),float ) ; U_L_x = np.zeros( (nx,M),float )
	U_0_xk = np.zeros( (nx,M),float ) ; U_L_xk = np.zeros( (nx,M),float )
	
	# singular end point locations and periodicity length for even extension
	a = x[0] ; b=x[-1]; P = 2*(b-a)    
	
	for i in range(nx):
		for j in range(M):
			n = 2*j+1						
			U_0[i,j] = U_n(x[i],a,P,n)                                  # function left
			U_L[i,j] = U_n(x[i],b,P,n)                                  # function right
			
			U_0_x[i,j] = kth_deriv_U_n(x[i],a,P,n,1)                    # 1st deriv left
			U_L_x[i,j] = kth_deriv_U_n(x[i],b,P,n,1)                    # 1st deriv right
			#U_0_x[i,j] = ddx_U_n(x[i],a,P,n)                    # 1st deriv left
			#U_L_x[i,j] = ddx_U_n(x[i],b,P,n)                    # 1st deriv right
			
			U_0_xk[i,j] = kth_deriv_U_n(x[i],a,P,n,k)                   # kth deriv left
			U_L_xk[i,j] = kth_deriv_U_n(x[i],b,P,n,k)                   # kth deriv right
			basis_functions=[U_0,U_L,U_0_x,U_L_x,U_0_xk,U_L_xk]
	return basis_functions


def ddx(f,x,frac,Q,LU,basis_functions):
#-------------------------------------------------------------------------------------------------
#    compute df/dx with f given at equally spaced grid points in [0,L] using the
#    Bernoulli-cosine approach to deal with the discontinuous derivatives at x=0,L
#    for the even extension of f into [0,2L)
#
#    frac is the fraction of the high wavenumber range to filter in cos transform differentiation
#         frac<0 ==> 2/3 rule, frac0 ==> no filtering
#    Q is the highest odd order of the U_n terms in the singular series expansions
#         (Q+1)/2 is the number of terms in the expansions
#    LU = [lu_A,lu_B,piv_A,piv_B] the prefactored matrices constraining the expansion coeffs
#
#    if Q=0, just do standard (filtered) term x term cosine differentiation
#-------------------------------------------------------------------------------------------------
	import numpy as np
	from Bernoulli_polynomials import solve_for_expansion_coeffs,U_series_expansion
	from Fourier_stuff import dst_filtered
	
	L = x[-1]-x[0]
	order = 1          # 1st derivative
	flag = -1          # even extend data, use cos series
	
	if( Q==0 ):
	
		dfdx = dst_filtered(f,order,L,flag,frac)
		
	else:
	
		s0=np.zeros_like(x); sL=np.zeros_like(x)
		s0_x=np.zeros_like(x); sL_x=np.zeros_like(x)
	
		M = (Q+1)/2
	
		# unpack  LU 
		[lu_A,lu_B,piv_A,piv_B] = LU
	
		# unpack  basis functions
		[U_0,U_L,U_0_x,U_L_x] = basis_functions
	
		#---------------------------------------------------------------
		#  solve the matrix problems for coeffs for 2 series expansions
		#---------------------------------------------------------------
		x0 = 0.
		A = solve_for_expansion_coeffs(lu_A,piv_A,f,Q,x0)
		x0 = L
		B = solve_for_expansion_coeffs(lu_B,piv_B,f,Q,x0)
	
		#---------------------------------------------------------------
		#  construct the 2 series for x in [0,L], add them to get f_s(x)
		#  also differentiate this series
		#---------------------------------------------------------------
		for j in range(M):
			s0   = s0   + A[j]*U_0[:,j]
			s0_x = s0_x + A[j]*U_0_x[:,j]
		
			sL   = sL   + B[j]*U_L[:,j]
			sL_x = sL_x + B[j]*U_L_x[:,j]
	 
		# combine the A and B series
		f_s   = s0   + sL 
		f_s_x = s0_x + sL_x  

		# extract the smooth, even extendable part that is well approximated by a cosine series
		f_Q = f - f_s    # should be Q times differentiable when even expanded 

		#---------------------------------------------------------------
		#  use standard cosine transform to differentiate f_Q
		#---------------------------------------------------------------
		f_Q_x = dst_filtered(f_Q,order,L,flag,frac)
	
		dfdx = f_Q_x + f_s_x
		
	return dfdx


def grad(f,x,z,frac,LU_x,LU_z,Q,basis_functions):
#-------------------------------------------------------------------------------------------------
#    compute grad f[x,z] = f_x and f_z using Bernoulli/Cosine scheme
#    f[x,z] w/ no particular symmetry on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(f)	
	f_x = np.zeros_like(f) ; f_z = np.zeros_like(f)
	[x_basis_functions,z_basis_functions] = basis_functions
	
	# compute f_x at each z
	for k in range(nz):
		f_x[:,k] = ddx(f[:,k],x,frac,Q,LU_x,x_basis_functions)
	
	# compute f_z at each x
	for i in range(nx):
		f_z[i,:] = ddx(f[i,:],z,frac,Q,LU_z,z_basis_functions)					
	return f_x,f_z

	
def divergence(u,w,x,z,frac,LU_x,LU_z,Q,basis_functions):
#-------------------------------------------------------------------------------------------------
#    compute div = u_x + w_z using Bernoulli/Cosine scheme
#    u[x,z] and w[x,z] w/ no particular symmetry on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(u)
	u_x = np.zeros_like(u) ; w_z = np.zeros_like(w)
	[x_basis_functions,z_basis_functions] = basis_functions
	
	# compute u_x
	for k in range(nz):
		u_x[:,k] = ddx(u[:,k],x,frac,Q,LU_x,x_basis_functions)
	
	# compute f_z
	for i in range(nx):
		w_z[i,:] = ddx(w[i,:],z,frac,Q,LU_z,z_basis_functions)
	
	div = u_x + w_z	
	return div


def grad2(f,x,z,frac,LU_x,LU_z,Q,basis_functions):
#-------------------------------------------------------------------------------------------------
#    compute grad^2 f[x,z] as div (grad f) using Bernoulli/Cosine scheme
#    f[x,z] w/ no particular symmetry on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	f_x,f_z = grad(f,x,z,frac,LU_x,LU_z,Q,basis_functions)
	grad2_f = divergence(f_x,f_z,x,z,frac,LU_x,LU_z,Q,basis_functions)
	return grad2_f


def endpoint_deriv(f,x,ii):
#-------------------------------------------------------------------------------------------------
#    use one sided FD scheme to compute the d/dx at an endpoint
#    if ii=0, estimate deriv at x[0], if ii=-1 estimate deriv at x[-1]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	dx = x[1]-x[0]
	nx = x.size
	
	if(ii == 0):
		i = 0
		df = ( -25.*f[i] + 48.*f[i+1] - 36.*f[i+2] + 16.*f[i+3] - 3.*f[i+4] )/(12.*dx)	
	elif(ii == -1):
		i = nx-1
		df = (  25.*f[i] - 48.*f[i-1] + 36.*f[i-2] - 16.*f[i-3] + 3.*f[i-4] )/(12.*dx)
	else:
		print("endpoint_deriv   ii must be 0 or -1 to indicate which end"); exit()
	return df


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

def exp(x,gamma):
#-----------------------------------------------------------------------------------------
#  Define functions that exponentially decay near endpoints of x array
#  the step functions have NONzero derivative at the ends and decay from 1 at the ends
#  to 0 in the interior of the domain over the scale gamma 
#-----------------------------------------------------------------------------------------
	import numpy as np
	a = x[0] ; b = x[-1]
	left  = np.exp(-((x-a)/gamma)) 
	right = np.exp(-((b-x)/gamma))
	return left,right

def create_2d_variables(nx,nz):
	import numpy as np
	u = np.zeros((nx,nz),dtype=float)
	v = np.zeros((nx,nz),dtype=float)
	w = np.zeros((nx,nz),dtype=float)
	b = np.zeros((nx,nz),dtype=float)
	eta = np.zeros((nx,nz),dtype=float)
	zeta = np.zeros((nx,nz),dtype=float)
	ustar = np.zeros((nx,nz),dtype=float)
	vstar = np.zeros((nx,nz),dtype=float)
	wstar = np.zeros((nx,nz),dtype=float)
	g2u = np.zeros((nx,nz),dtype=float)	
	g2v = np.zeros((nx,nz),dtype=float)
	g2w = np.zeros((nx,nz),dtype=float)
	g2b = np.zeros((nx,nz),dtype=float)
	vars = [u,v,w,b,eta,zeta,ustar,vstar,wstar,g2u,g2v,g2w,g2b]
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


def construct_psi_u(x,z,BVALS):
	import numpy as np		
	nx = x.size ; nz=z.size
	psi_u = np.zeros((nx,nz),float)
	# unpack the boundary values into 1d arrays and get ustar and wstar
	[U_east,U_west,U_z_bot,U_z_top,ustar,wstar] = BVALS
		
	#------------------------------------------------------------------------------	
	# construct an array that matches the required Dirichlet bcs at east/west
	#------------------------------------------------------------------------------
	npts = 3        # size of x transition layer at E/W
	gamma = npts*( x[1]-x[0] )
	step_e,step_w = step(x,gamma)	
	for k in range(nz):
		a = (U_east[k]-ustar[ 0,k])
		b = (U_west[k]-ustar[-1,k])	
		psi_u[:,k] = a*step_e[:] + b*step_w[:]
	return psi_u


def construct_psi_w(x,z,BVALS):
	import numpy as np
	nx = x.size ; nz=z.size
	psi_w = np.zeros((nx,nz),float)
	# unpack the boundary values into 1d arrays and get ustar and wstar
	[W_bot,W_top,W_x_east,W_x_west,ustar,wstar] = BVALS
		
	#------------------------------------------------------------------------------	
	# construct an array that matches the required Dirichlet bcs at top/bottom
	#------------------------------------------------------------------------------
	npts = 6     # size of z transition layer at top/bottom
	gamma = npts*( z[1]-z[0] )
	step_b,step_t = step(z,gamma)	
	for i in range(nx):
		a = (W_bot[i]-wstar[i, 0])
		b = (W_top[i]-wstar[i,-1])	
		psi_w[i,:] = a*step_b[:] + b*step_t[:]
	return psi_w
	
	
def apply_bcs(u,v,w,b,x,z,BVALS):
	import numpy as np
	nx = x.size  ; nz = z.size
	
	# unpack boundary information into 1d arrays
	[BVALS_EW,BVALS_BT] = BVALS
	[U_east,U_west,V_east,V_west,W_east,W_west,B_east,B_west,W_x_east,W_x_west,W_xx_east,W_xx_west] = BVALS_EW
	[U_bot,U_top,V_bot,V_top,W_bot,W_top,B_bot,B_top,U_z_bot,U_z_top,U_zz_bot,U_zz_top,W_z_bot,W_z_top,W_zz_bot,W_zz_top] = BVALS_BT
	
	# correct w in BL near E/W boundaries; ; match functions at x[ii]
	ii=2*3
	for k in range(nz):
		f = W_east[k] + W_x_east[k]*x[:] #+ (1./2.)*W_xx_east[k]*x[:]**2                  # low order Taylor series expansion
		offset = w[ii-1,k]-f[ii-1]
		w[0:ii,k] = f[0:ii] + offset
		
		f = W_west[k] + W_x_west[k]*(x[:]-x[-1]) #+ (1./2.)*W_xx_west[k]*(x[:]-x[-1])**2  # low order Taylor series expansion
		offset = w[nx-ii-1,k]-f[nx-ii-1]
		w[-ii:nx,k] = f[-ii:nx] + offset
	
	# correct u in BL near B/T boundaries; match functions at z[kk]
	kk=6
	for i in range(nx):
		f = U_bot[i] + U_z_bot[i]*z[:] #+ (1./2.)*U_zz_bot[i]*z[:]**2
		offset = u[i,kk-1] - f[kk-1]
		u[i,0:kk] = f[0:kk] #+ offset
		
		f = U_top[i] + U_z_top[i]*(z[:]-z[-1]) #+ (1./2.)*U_zz_top[i]*(z[:]-z[-1])**2
		offset = u[i,nz-kk-1]-f[nz-kk-1]
		u[i,-kk:nz] = f[-kk:nz] #+ offset
		
		f = W_bot[i] + W_z_bot[i]*z[:] #+ (1./2.)*W_zz_bot[i]*z[:]**2
		offset = w[i,kk-1] - f[kk-1]
		w[i,0:kk] = f[0:kk] #+ offset
		
		f = W_top[i] + W_z_top[i]*(z[:]-z[-1]) #+ (1./2.)*W_zz_top[i]*(z[:]-z[-1])**2
		offset = w[i,nz-kk-1]-f[nz-kk-1]
		w[i,-kk:nz] = f[-kk:nz] #+ offset
	
	# insert the scalar boundary values
	for k in range(nz):
		v[0,k] = V_east[k] ; v[-1,k] = V_west[k]
		b[0,k] = B_east[k] ; b[-1,k] = B_west[k]
	
	for i in range(nx):	
		v[i,0] = V_bot[i]  ; v[i,-1] = V_top[i]	
		b[i,0] = B_bot[i]  ; b[i,-1] = B_top[i]	
	
	#  smooth over scale npts*dz near the bottom and top boundaries
	dir = 'z' ; npts = 10
	for i in range(3):
		u = boundary_smooth(u,x,z,dir,npts)
		w = boundary_smooth(w,x,z,dir,npts)
		b = boundary_smooth(b,x,z,dir,npts)
		
	return u,v,w,b	


def boundary_smooth(f,x,z,dir,npts):
	import numpy as np
	nx = x.size  ;  nz = z.size
	dx = x[1]-x[0] ; dz = z[1]-z[0]
	
	if( dir=='x' or dir=='both' ):
		gamma = npts*dx
		for k in range(nz):
			f[:,k] = near_boundary_smoothing(f[:,k],x,gamma)
	
	if( dir=='z' or dir=='both' ):
		gamma = npts*dz
		for i in range(nx):
			f[i,:] = near_boundary_smoothing(f[i,:],z,gamma)
	return f
	

def near_boundary_smoothing(f,x,gamma):
	import numpy as np
	nx = x.size ; fs = f.copy()
	w = np.exp(-x/gamma) + np.exp(-(x[-1]-x)/gamma)  # 1 at boundary --> 0 in interior
	for i in range(1,nx-1):
		fs[i] = (w[i]*f[i-1] + f[i] + w[i]*f[i+1]) / (2*w[i]+1.)		
	return fs	

def add_body_force(sigma,t,x,z,A,rhs):
#  add an oscillating point source forcing term to rhs
	import numpy as np
	nx=x.size                       ;   nz=z.size
	dx = x[1]-x[0] ; gamma_x=2*dx   ;   dz = z[1]-z[0] ; gamma_z=2*dz  
	x0 = (x[-1]-x[0])/2.            ;   z0 = (z[-1]-z[0])/2.
	xx = A*np.cos(sigma*t)   # amplitude and time dependence
	for k in range(nz):
		for i in range(nx):
			W = np.exp(-((x[i]-x0)/gamma_x)**2) * np.exp(-((z[k]-z0)/gamma_z)**2)  # spatial shape
			rhs[i,k] = rhs[i,k] + W*xx
	return rhs




def time_step(n,dt,U,RHS,AB_order,idx):
#  n is the current time step, time stepping starts at n=0
#  time step the eqns one step using Adams Bashforth schemes
#  toggle the indices "idx"
#  RHS[nx,nz,var_id,4]     the last index is for the 4 time levels required for AB4
	import numpy as np
	
	inc = 5
	[u,v,w,b,eta,zeta] = U
	[MM0,MM1,MM2,MM3,method] = idx
	if( n < AB_order or np.mod(n,inc)==0 ): print("Step  ",n,method)
	
	
	if(method=='euler'):      # Euler
		u = u + dt*RHS[:,:,0,MM0]
		v = v + dt*RHS[:,:,1,MM0]
		w = w + dt*RHS[:,:,2,MM0]
		b = b + dt*RHS[:,:,3,MM0]
		eta = eta + dt*RHS[:,:,4,MM0]
		zeta = zeta + dt*RHS[:,:,5,MM0]
		if( AB_order > 1 ): method = 'AB2'
                        	
	if(n > 0 and method=='AB2'):       # AB2
		dt2 = dt/2.
		u = u + dt2*( 3.*RHS[:,:,0,MM0] - RHS[:,:,0,MM1] )
		v = v + dt2*( 3.*RHS[:,:,1,MM0] - RHS[:,:,1,MM1] )
		w = w + dt2*( 3.*RHS[:,:,2,MM0] - RHS[:,:,2,MM1] )
		b = b + dt2*( 3.*RHS[:,:,3,MM0] - RHS[:,:,3,MM1] )
		eta = eta + dt2*( 3.*RHS[:,:,4,MM0] - RHS[:,:,4,MM1] )
		zeta = zeta + dt2*( 3.*RHS[:,:,5,MM0] - RHS[:,:,5,MM1] )
		if( AB_order > 2 ): method = 'AB3'
		
	if(n > 1 and method=='AB3'):       # AB3
		dt12 = dt/12.
		u = u + dt12*( 23.*RHS[:,:,0,MM0] - 16.*RHS[:,:,0,MM1] + 5.*RHS[:,:,0,MM2] )
		v = v + dt12*( 23.*RHS[:,:,1,MM0] - 16.*RHS[:,:,1,MM1] + 5.*RHS[:,:,1,MM2] )
		w = w + dt12*( 23.*RHS[:,:,2,MM0] - 16.*RHS[:,:,2,MM1] + 5.*RHS[:,:,2,MM2] )
		b = b + dt12*( 23.*RHS[:,:,3,MM0] - 16.*RHS[:,:,3,MM1] + 5.*RHS[:,:,3,MM2] )
		eta = eta + dt12*( 23.*RHS[:,:,4,MM0] - 16.*RHS[:,:,4,MM1] + 5.*RHS[:,:,4,MM2] )
		zeta = zeta + dt12*( 23.*RHS[:,:,5,MM0] - 16.*RHS[:,:,5,MM1] + 5.*RHS[:,:,5,MM2] )
		if( AB_order > 3 ): method = 'AB4'
		
	if(n > 2 and method=='AB4'):        # AB4  
		dt24 = dt/24.
		u = u + dt24*( 55.*RHS[:,:,0,MM0] - 59.*RHS[:,:,0,MM1] + 37.*RHS[:,:,0,MM2] - 9.*RHS[:,:,0,MM3] )
		v = v + dt24*( 55.*RHS[:,:,1,MM0] - 59.*RHS[:,:,1,MM1] + 37.*RHS[:,:,1,MM2] - 9.*RHS[:,:,1,MM3] )
		w = w + dt24*( 55.*RHS[:,:,2,MM0] - 59.*RHS[:,:,2,MM1] + 37.*RHS[:,:,2,MM2] - 9.*RHS[:,:,2,MM3] )
		b = b + dt24*( 55.*RHS[:,:,3,MM0] - 59.*RHS[:,:,3,MM1] + 37.*RHS[:,:,3,MM2] - 9.*RHS[:,:,3,MM3] )
		eta = eta + dt24*( 55.*RHS[:,:,4,MM0] - 59.*RHS[:,:,4,MM1] + 37.*RHS[:,:,4,MM2] - 9.*RHS[:,:,4,MM3] )
		zeta = zeta + dt24*( 55.*RHS[:,:,5,MM0] - 59.*RHS[:,:,5,MM1] + 37.*RHS[:,:,5,MM2] - 9.*RHS[:,:,5,MM3] )
	
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
 
	U = [u,v,w,b,eta,zeta]                     # updated dependent variables
 	jdx = [MM0,MM1,MM2,MM3,method]    # updated index array
	return U,jdx
	
def high_order_ddx(f,x,frac,Q,LU,basis_functions,k,inc):
#-------------------------------------------------------------------------------------------------
#    compute df/dx with f given at equally spaced grid points in [0,L] using the
#    Bernoulli-cosine approach to deal with the discontinuous derivatives at x=0,L
#    for the even extension of f into [0,2L)
#
#    k is order of the derivative order, assume basis functions available 
#
#    frac is the fraction of the high wavenumber range to filter in cos transform differentiation
#         frac<0 ==> 2/3 rule, frac0 ==> no filtering
#    Q is the highest odd order of the U_n terms in the singular series expansions
#         (Q+1)/2 is the number of terms in the expansions
#    LU = [lu_A,lu_B,piv_A,piv_B] the prefactored matrices constraining the expansion coeffs
#
#    if Q=0, just do standard (filtered) term x term cosine differentiation
#-------------------------------------------------------------------------------------------------
	import numpy as np
	import matplotlib
	import matplotlib.pyplot as plt
	from Bernoulli_polynomials import solve_for_expansion_coeffs,U_series_expansion
	from Fourier_stuff import dst_filtered
	
	L = x[-1]-x[0]
	flag = -1          # even extend data, use cos series
	nx = x.size
	dx = x[1]-x[0]
	
	if( Q==0 ):
	
		dfdx = dst_filtered(f,k,L,flag,frac)
		
	else:
	
		s0=np.zeros_like(x); sL=np.zeros_like(x)
		s0_x=np.zeros_like(x); sL_x=np.zeros_like(x)
		s0_xk=np.zeros_like(x); sL_xk=np.zeros_like(x)
	
		M = (Q+1)/2
	
		# unpack  LU 
		[lu_A,lu_B,piv_A,piv_B] = LU
	
		# unpack  basis functions
		[U_0,U_L,U_0_x,U_L_x,U_0_xk,U_L_xk] = basis_functions
	
		#---------------------------------------------------------------
		#  solve the matrix problems for coeffs for 2 series expansions
		#---------------------------------------------------------------
		x0 = 0.
		A = solve_for_expansion_coeffs(lu_A,piv_A,f,Q,x0,inc)
		x0 = L
		B = solve_for_expansion_coeffs(lu_B,piv_B,f,Q,x0,inc)
	
		#---------------------------------------------------------------
		#  construct the 2 series for x in [0,L], add them to get f_s(x)
		#  also differentiate this series k times
		#---------------------------------------------------------------
		for j in range(M):
			s0    = s0    + A[j]*U_0[:,j]
			s0_x  = s0_x  + A[j]*U_0_x[:,j]
			s0_xk = s0_xk + A[j]*U_0_xk[:,j]
		
			sL    = sL    + B[j]*U_L[:,j]
			sL_x  = sL_x  + B[j]*U_L_x[:,j]
			sL_xk = sL_xk + B[j]*U_L_xk[:,j]
			
		# combine the A and B series
		f_s    = s0    + sL             # left and right series matching behavior of f near ends
		f_s_x  = s0_x  + sL_x           # 1st derivative of these series
		f_s_xk = s0_xk + sL_xk          # kth derivative of these series

		# extract the smooth, even extendable part that is well approximated by a cosine series
		f_Q = f - f_s    # should be ~Q times differentiable when even expanded  (it's not)
		
		#plt.plot(x,f,'k',x,f_Q,'r',x,f_s,'b')
		plt.plot(x,f_s_xk)
		plt.show()
		

		#---------------------------------------------------------------
		#  use standard cosine transform to differentiate f_Q
		#---------------------------------------------------------------
		f_Q_x = dst_filtered(f_Q,k,L,flag,frac)
		
		if( k==1 ):
			kth_deriv =  f_Q_x + f_s_x
		elif( k>1 ):
			kth_deriv = f_Q_x + f_s_xk
			#gamma = L/10 ; delta = dx*4
			#for i in range(6):
			#	kth_deriv = near_boundary_smoothing2(kth_deriv,x,gamma)
			#	kth_deriv = to_zero_at_ends(kth_deriv,x,delta)
				
	return kth_deriv

def near_boundary_smoothing2(f,x,gamma):
	import numpy as np
	nx = x.size ; fs = f.copy()
	
	# avg end values 
	fs[0] = (f[0]+f[1])/2.
	fs[nx-1] = (f[nx-1]+f[nx-2])/2.
	
	w = np.exp(-x/gamma) + np.exp(-(x[-1]-x)/gamma)  # 1 at boundary --> 0 in interior
	for i in range(1,nx-1):
		fs[i] = (w[i]*f[i-1] + f[i] + w[i]*f[i+1]) / (2*w[i]+1.)
			
	return fs	

	
def shift_smooth_unshift(f,x,gamma):
	import numpy as np
	nx = x.size ; fs = f.copy()
	p = 2
	
	# left side
	fs = f[:] - f[0]
	for i in range(nx):
		w = 1. - np.exp(-((x[i]-x[0])/gamma)**p)
		fs[i]=f[i]*w
	f[:] = fs[:] + f[0]
	
	# right side
	fs = f[:] - f[-1]
	for i in range(nx):
		w = 1. - np.exp(-((x[i]-x[-1])/gamma)**p)
		fs[i]=f[i]*w
	f[:] = fs[:] + f[-1]			
	return fs	

def to_zero_at_ends(f,x,gamma):
	import numpy as np
	nx = x.size 
	p = 4	
	# left side
	for i in range(nx):
		w = 1. - np.exp(-((x[i]-x[0])/gamma)**p) - np.exp(-((x[i]-x[-1])/gamma)**p)
		f[i]=f[i]*w	
	return f	

