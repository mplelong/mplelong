#-----------------------------------------------------------------------------------------
#  KBW python scripts for stuff related to Bernoulli polynomials
#    
#    -----------------------------
#      list of functions defined:
#    -----------------------------
#-----------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
	
	
def even_extend(f):
#---------------------------------------------------------------------------------
#    f(x) at N uniformly spaced discrete grid points in closed interval [0,L] 
#    f_even  the even extension of f into the open interval [0,2L], (M, periodic)
#---------------------------------------------------------------------------------
	import numpy as np
	N = f.size; M = (N-1)*2
	f_even = np.concatenate([f,  f[1:-1][::-1]])   # even extension of data
	return f_even

	
def odd_extend(f):
#---------------------------------------------------------------------------------
#    f(x) at N uniformly spaced discrete grid points in closed interval [0,L] 
#    f_odd  the odd extension of f into the open interval [0,2L], (M, periodic)
#---------------------------------------------------------------------------------
	import numpy as np
	N = f.size; M = (N-1)*2
	f_odd = np.concatenate([f,  -f[1:-1][::-1]])   # odd extension of data
	return f_odd

	
def B(n):
#--------------------------------------------------------------------------------------
#    Evaluate the nth Bernoulli number using the summation formula (2) in
#    Komatsu and Pita Ruiz V. Math Commun. 21, 2016 pp. 127-140
#--------------------------------------------------------------------------------------
	import numpy as np
	from scipy.special import comb  #  comb(n,k)  "n choose k" = n!/( k!(n-k)! )
	sum = 0.
	for i in range(n+1):
		for j in range(i+1):
			sum = sum + (-1)**(i+j) * comb(n+1,j)/comb(n,i) * (i-j)**n
	sum = sum/(n+1.)					
	return sum

	
def B_n(x,n):
#--------------------------------------------------------------------------------------
#    Evaluate the nth Bernoulli polynomial B_n(x) using the 
#    "explicit formula" in Wikepedia https://en.wikipedia.org/wiki/Bernoulli_polynomials	  
#    ==> much more accurate than Komatsu and Pita Ruiz V. Math Commun. 21, 2016 pp. 127-140
#--------------------------------------------------------------------------------------
	import numpy as np
	from scipy.special import comb  #  comb(n,k)  "n choose k" = n!/( k!(n-k)! )
	
	sum = 0.	
	for k in range(n+1):
		sum = sum + comb(n,k)*B(n-k)*x**k				
	return sum


def ddx_B_n(x,n):
#-------------------------------------------------------------------------------
#    Evaluate d/dx of the nth Bernoulli polynomial B_n(x) 
#-------------------------------------------------------------------------------
	import numpy as np
	ddx = n*B_n(x,n-1)    # d/dx B_n(x) = n*B_n-1(x)				
	return ddx


def kth_deriv_B_n(x,n,k):
#-------------------------------------------------------------------------------
#    Evaluate kth derivative of the nth Bernoulli polynomial B_n(x)
#    input params: n is Bernoulli polynomial index/order 
#                  m is the order of the derivative to evaluate at position x  
#-------------------------------------------------------------------------------
	import numpy as np
	from scipy.special import factorial as fac     #   N.B. fac(integer) returns a float
	if( n >= k ):
		deriv = fac(n)/fac(n-k) * B_n(x,n-k)   #  last term is B_(n-k)(x)
	else:
		deriv = 0.	
	return deriv

def U_n(x,a,P,n):
#--------------------------------------------------------------------------------------
#    Evaluate the function U_n(x-a) defined in Eckhoff 98
#    MATHEMATICS OF COMPUTATION
#    Volume 67, Number 223, July 1998, pp. 1063-1087
#    The periodic in [0,1) function U_n(x) is a piecewise polynomial of degree n + 1
#     
#    a is shift, P is periodicity length for U_n(x), n is the order
#    N.B.  argument of the Bernoulli polynomials must be between 0 and 1
#--------------------------------------------------------------------------------------
	import numpy as np
	from scipy.special import factorial as fac
	xsi = (x-a)/P
	#c = -(P)**n / fac(n+1)
	c = -1.
	#  wrap xsi as necessary for periodicity
	if( xsi < 0 ): 
		xsi = 1. + xsi	
	elif( xsi > 1 ):
		xsi = xsi - 1.
	# endpoint should have xsi=1 not xsi = 0 
	if( x==P/2. and a==P/2. ): xsi = 1.0	
	ans = c*B_n(xsi,n+1)					
	return ans
	
def ddx_U_n(x,a,P,n):
#--------------------------------------------------------------------------------------
#    use the differentiation formula for the Bernoulli polynomials
#
#    Evaluate the function U_n(x-a) defined in Eckhoff 98
#    MATHEMATICS OF COMPUTATION
#    Volume 67, Number 223, July 1998, pp. 1063-1087
#    The periodic in [0,1) function U_n(x) is a piecewise polynomial of degree n + 1
#     
#    a is shift, L is periodicity length for U_n(x), n is the order
#    N.B.  argument of the Bernoulli polynomials must be between 0 and 1
#
#  U_n = c*B_n(xsi,n+1)                        U_n in terms of B_n+1	, xsi = (x-a)/L
#  d/dx U_n = c * d/dxsi B_n+1(xsi) * dxsi/dx
#  d/dxsi B_n+1(xsi) = (n+1)*B_n(xsi)
#--------------------------------------------------------------------------------------
	import numpy as np
	from scipy.special import factorial as fac
	xsi = (x-a)/P
	dxsi_dx = 1./P
	#c = -(P)**n / fac(n+1)
	c = -1.
	
	#  wrap xsi as necessary for periodicity
	if( xsi < 0 ): 
		xsi = 1. + xsi	
	elif( xsi > 1 ):
		xsi = xsi - 1.			
	# endpoint should have xsi=1 not xsi = 0 
	if( x==P/2. and a==P/2. ): xsi = 1.0 
				
	ddxsi = c * (n+1.)*B_n(xsi,n)    # d/dxsi of U_n(xsi)
	ddx = ddxsi * dxsi_dx            # d/dx of U_n(x)
	return ddx
	
def kth_deriv_U_n(x,a,P,n,k):
#--------------------------------------------------------------------------------------
#    use the differentiation formula for the Bernoulli polynomials
#
#    Evaluate the function U_n(x-a) defined in Eckhoff 98
#    MATHEMATICS OF COMPUTATION
#    Volume 67, Number 223, July 1998, pp. 1063-1087
#    The periodic in [0,1) function U_n(x) is a piecewise polynomial of degree n + 1
#     
#    a is shift, L is periodicity length for U_n(x), n is the order
#    N.B.  argument of the Bernoulli polynomials must be between 0 and 1
#
#  U_n = c*B_n(xsi,n+1)                        U_n in terms of B_n+1	, xsi = (x-a)/L
#  d/dx U_n = c * d/dxsi B_n+1(xsi) * dxsi/dx
#  d/dxsi B_n+1(xsi) = (n+1)*B_n(xsi)            ===> apply k times
#--------------------------------------------------------------------------------------
	import numpy as np
	from scipy.special import factorial as fac
	xsi = (x-a)/P
	dxsi_dx = 1./P
	#c = -(P)**n / fac(n+1)
	c = -1.
	
	#  wrap xsi as necessary for periodicity
	if( xsi < 0 ): 
		xsi = 1. + xsi	
	elif( xsi > 1 ):
		xsi = xsi - 1.			
	# endpoint should have xsi=1 not xsi = 0 
	if( x==P/2. and a==P/2. ): xsi = 1.0 
	
	if( n >= k ):
		ddxsi = c * fac(n)/fac(n-k) * B_n(xsi,n-k)  # kth deriv of U_n(xsi) wrt to xsi
		deriv = ddxsi * dxsi_dx**k                  # kth deriv of U_n(x) wrt x
	else:
		deriv = 0.			
	return deriv
	
	
def setup_factor_matrix(x,Q,x0):
#--------------------------------------------------------------------------------------
#    set up the matrix problem that determines the weights of the U_n(x) functions
#    x is a vector of equally spaced points [0,L] where f is defined
#    the even extension of f is periodic over length 2L. 
#    x0 is the endpoint around which the U_n functions are expanded 
#    (i.e. x=0 and x=L) for the A and B series expansions. The matrix eqns
#    constrain the Q term expansions to match f at  (Q+1)/2 points near 0 and L.
#--------------------------------------------------------------------------------------
	import numpy as np
	from scipy.linalg import lu_factor
	
	nx = x.size
	npts = (Q+1)/2  ;  M = np.zeros((npts,npts))
	
	if(x0==x[0]):
		istart = 0                  # array index of first x val to include
	elif(x0==x[-1]): 
		istart = nx - npts         # array index of first x val to include
	else:
		print("setup_factor_matrices: problem with coordinate shift x0",x0)
		exit()
			
	a = x0            # shift for evaluating U_n(x-a)    either 0 or L
	L = x[-1]-x[0]    # original domain length
	P = 2*L           # periodicity length for even extensions
	
	# fill the elements of the matrix M
	for i in range(npts):                 # matrix rows, x location varies
		xval = x[istart+i]
		for j in range(npts):             # matrix cols, odd polynomial order varies
			n = 2*j+1
			M[i,j] = U_n(xval,a,P,n)
	
	# decompose M into LU		
	lu,piv = lu_factor(M)
	return lu,piv
	
def solve_for_expansion_coeffs(lu,piv,f,Q,x0):
#--------------------------------------------------------------------------------------
#    fill the rhs vector and solve for the expansion coefficients for the U_n series
#    if the expansion point x0=0, then the coeffs are A1,A3,..A_Q and the series
#    should match f at the (Q+1)/2 points near x=0
#    if the expansion point x0=L, then the coeffs are B1,B3,..B_Q and the series
#    should match f at the (Q+1)/2 points near x=L
#--------------------------------------------------------------------------------------
	import numpy as np
	from scipy.linalg import lu_solve
	
	M = (Q+1)/2  ; rhs = np.zeros((M,1))
	nx = f.size
	if( x0==0. ):
		istart = 0                 # array index of first x val to include
	else:
		istart = nx - M            # array index of first x val to include

	# fill the rhs vector
	for i in range(M):
		rhs[i] = f[istart+i]
		
	coeffs = lu_solve((lu, piv), rhs)
	return coeffs


def U_series_expansion(x,Q,x0,coeffs):
#--------------------------------------------------------------------------------------
#    x is an nx vector of equally spaced points in [0,L]
#    construct the U_n series expansion with basis functions singular at x=x0
#    given x, Q the highest order in the (Q+1)/2 odd term expansion and 
#    the expansion coefficients coeffs[:]. Also construct the first derivative.
#--------------------------------------------------------------------------------------
	import numpy as np
	
	L = x[-1] - x[0]     # domain length
	P = 2*L              # periodicity scale for even extensions
	nx = x.size
	s = np.zeros_like(x) ; s_x = np.zeros_like(x)
	
	for i in range(nx):
		for j in range( (Q+1)/2 ):
			n = 2*j + 1
			s[i]   = s[i]   + coeffs[j]*U_n(x[i],x0,P,n)
			s_x[i] = s_x[i] + coeffs[j]*ddx_U_n(x[i],x0,P,n)
	return s,s_x


def U_series_expansion_k(x,Q,x0,coeffs,k):
#--------------------------------------------------------------------------------------
#    x is an nx vector of equally spaced points in [0,L]
#    construct the U_n series expansion with basis functions singular at x=x0
#    given x, Q the highest order in the (Q+1)/2 odd term expansion and 
#    the expansion coefficients coeffs[:]. 
#    Also construct the first derivative and kth derivatives.
#--------------------------------------------------------------------------------------
	import numpy as np
	
	L = x[-1] - x[0]     # domain length
	P = 2*L              # periodicity scale for even extensions
	nx = x.size
	s = np.zeros_like(x) ; s_x = np.zeros_like(x); s_xk = np.zeros_like(x)
	
	for i in range(nx):
		for j in range( (Q+1)/2 ):
			n = 2*j + 1
			s[i]    = s[i]    + coeffs[j]*U_n(x[i],x0,P,n)
			s_x[i]  = s_x[i]  + coeffs[j]*ddx_U_n(x[i],x0,P,n)
			s_xk[i] = s_xk[i] + coeffs[j]*kth_deriv_U_n(x[i],x0,P,n)
	return s,s_x,s_xk















 
