def prepare_BernoulliCosine(x,Q):
#-------------------------------------------------------------------------------------------------
#    construct and factor the matrices required for the Bernoulli expansion coefficients
#    LU = [lu_A,lu_B,piv_A,piv_B] the prefactored matrices constraining the expansion coeffs
#-------------------------------------------------------------------------------------------------
	import numpy as np
	from Bernoulli_polynomials import setup_factor_matrix
	x0 = x[0]     # even extension of f has discontinuous deriv at x=0
	lu_A,piv_A = setup_factor_matrix(x,Q,x0)
	x0 = x[-1]    # even extension of f also has discontinuous deriv at x=L
	lu_B,piv_B = setup_factor_matrix(x,Q,x0)
	LU = [lu_A,lu_B,piv_A,piv_B]
	return LU

def evaluate_basis_functions(x,Q):
#-------------------------------------------------------------------------------------------------
#    Construct the U_n basis functions for n=1,3,5... and their derivatives
#    at each x point for the (Q+1)/2 term series
#    The S_0 series are expended about x=0 while the S_L series are expanded about x=x[-1]=L
#-------------------------------------------------------------------------------------------------
	import numpy as np
	from Bernoulli_polynomials import U_n,ddx_U_n
	nx = x.size
	# number of terms in the series expansions
	M = (Q+1)/2
	U_0 = np.zeros( (nx,M),float )   ; U_L = np.zeros( (nx,M),float )
	U_0_x = np.zeros( (nx,M),float ) ; U_L_x = np.zeros( (nx,M),float )
	
	# singular end point locations and periodicity length for even extension
	a = x[0] ; b=x[-1]; P = 2*(b-a)    
	
	for i in range(nx):
		for j in range(M):
			n = 2*j+1						
			U_0[i,j] = U_n(x[i],a,P,n)
			U_L[i,j] = U_n(x[i],b,P,n)
			U_0_x[i,j] = ddx_U_n(x[i],a,P,n)
			U_L_x[i,j] = ddx_U_n(x[i],b,P,n)
			basis_functions=[U_0,U_L,U_0_x,U_L_x]
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


def grad_2d(f,x,z,frac,LU_x,LU_z,Q,basis_functions):
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

	
def divergence_2d(u,w,x,z,frac,LU_x,LU_z,Q,basis_functions):
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


def grad2_2d(f,x,z,frac,LU_x,LU_z,Q,basis_functions):
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

