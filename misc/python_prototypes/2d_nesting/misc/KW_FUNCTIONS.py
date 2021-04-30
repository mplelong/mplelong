#-----------------------------------------------------------------------------------------
#  KBW python scripts for general purpose use
#    
#    -----------------------------
#      list of functions defined:
#    -----------------------------
#     lu,piv,B = build_compact_matrices_1(n,L)                       build/factor compact differentiation matrices for 1st derivative
#     lu,piv,B = build_compact_matrices_2(n,L)                       build/factor compact differentiation matrices for 2nd derivative
#     df = compact_deriv(f,lu,piv,B)                                 compute f'[:] using compact schemes, direction/order set by input matrices
#     f_x,f_z = grad(f,lu_x,piv_x,B_x,lu_z,piv_z,B_z)                compute gradient of f[x,z] in [0,Lx],[0,Lz] 
#     div = divergence(u,w,lu_x,piv_x,B_x,lu_z,piv_z,B_z)            compute div = u_x + w_z    in [0,Lx],[0,Lz] 
#     grad2_f = grad2(f,lu_x,piv_x,B_x,lu_z,piv_z,B_z)               compute f_xx + f_zz        in [0,Lx],[0,Lz] 
#     inversion_matrix = create_inversion_matrix(nx,nz,Lx,Lz)        fill 2d array -1/(k^2 + m^2) for even-extended, periodic data
#     p = poisson_invert_2d(f,L)                                     solve p_xx + p_zz = f with homog. Neumann BCs
#     write_netcdf_2d(arrays,xvec,zvec,tval,filename,varnames,dims)
#
#     df = cos_deriv_filtered(f,n,L,frac)
#     high_even_deriv
#-----------------------------------------------------------------------------------------

	
	
	
	
def grad(f,lu_x,piv_x,B_x,lu_z,piv_z,B_z):
#-------------------------------------------------------------------------------------------------
#    compute grad f[x,z] = f_x and f_z using compact schemes
#    f[x,z] w/ no particular symmetry on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	#
	nx,nz = np.shape(f)	
	f_x = np.zeros_like(f) ; f_z = np.zeros_like(f)
	
	# compute f_x
	for k in range(nz):
		f_x[:,k] = compact_deriv(f[:,k],lu_x,piv_x,B_x)
	
	# compute f_z
	for i in range(nx):
		f_z[i,:] = compact_deriv(f[i,:],lu_z,piv_z,B_z)	
				
	return f_x,f_z

	
def divergence(u,w,lu_x,piv_x,B_x,lu_z,piv_z,B_z):
#-------------------------------------------------------------------------------------------------
#    compute div = u_x + w_z using compact schemes
#    u[x,z] and v[x,z] w/ no particular symmetry on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(u)
	u_x = np.zeros_like(u) ; w_z = np.zeros_like(w)
	
	# compute u_x
	for k in range(nz):
		u_x[:,k] = compact_deriv(u[:,k],lu_x,piv_x,B_x)
	
	# compute f_z
	for i in range(nx):
		w_z[i,:] = compact_deriv(w[i,:],lu_z,piv_z,B_z)
	
	div = u_x + w_z	
	return div


def grad2(f,lu_x,piv_x,B_x,lu_z,piv_z,B_z):
#-------------------------------------------------------------------------------------------------
#    compute grad^2 f[x,z] = f_xx + f_zz using compact schemes
#    input matrices constructed for 2nd derivatives
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(f)
	f_xx= np.zeros_like(f) ; f_zz = np.zeros_like(f)
	
	# compute f_xx
	for k in range(nz):
		f_xx[:,k] = compact_deriv(f[:,k],lu_x,piv_x,B_x)
		
	# compute f_zz
	for i in range(nx):
		f_zz[i,:] = compact_deriv(f[i,:],lu_z,piv_z,B_z)
		
	grad2_f = f_xx + f_zz
	
	return grad2_f





def create_inversion_matrix(nx,nz,Lx,Lz):
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
        
	MX = (nx-1)*2 ; MZ = (nz-1)*2    # size of even extended data array, Fourier wavnumber arrays
    
	dk = np.pi/Lx
	k = dk * np.array( range(MX/2+1) + range(-MX/2+1,0,1) )  # Fourier x wavenumbers for extended array
    
    
	dm = np.pi/Lz
	m = dm * np.array( range(MZ/2+1) + range(-MZ/2+1,0,1) )  # Fourier z wavenumbers for extended array

    
	inversion_matrix = np.zeros((MX,MZ),dtype=float)
	for kk in range(MZ):
		for ii in range(MX):
			denom = (k[ii]**2 + m[kk]**2)
			if( denom==0. ):
				inversion_matrix[ii,kk] = 0.0
			else:
				inversion_matrix[ii,kk] = -1./denom
    			
	return  inversion_matrix                           # return the 2d inversion matrix for later, repeated use


def poisson_invert_2d(f,inversion_matrix):
#-----------------------------------------------------------------------------------------
#  solve p_xx + p_zz = f(x,z)  
#      +  homogeneous Neumann conditions on the boundary of [0,Lx] x [0,Lz]
#  p and f are cosine expandable functions in both x and z directions 
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft2, ifft2
        
    NX,NZ = np.shape(f)                                      # size of input array
    MX = (NX-1)*2 ; MZ = (NZ-1)*2                            # size of even extended array        
    F = np.zeros((MX,MZ))                                    # array for holding 2d even extended data
    
    #-------------------------------------------------------
    # even extend data in x
    #-------------------------------------------------------
    for kk in range(NZ):
    	tmp = f[:,kk]
    	F[:,kk] = np.concatenate([tmp,  tmp[1:-1][::-1]])
    
    #-------------------------------------------------------
    # even extend data in z
    #-------------------------------------------------------
    for ii in range(MX):
    	tmp = F[ii,0:NZ]
    	F[ii,:] = np.concatenate([tmp,  tmp[1:-1][::-1]])
    
    
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
    

def write_netcdf_2d(arrays,xvec,zvec,tval,filename,varnames,dims):
#
#	Simple routine to write 2d spatial data to the netcdf file filename: 
#
#     arrays is a list of 2d data arrays
#     2d data stored in python arrays with layout f[x,z] with names varnames and dimensions dims
#     The 1d arrays x,z are assumed to be spatial arrays in [m]
#
#
#   5 flavors of netcdf files:
#     (NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET, NETCDF3_64BIT_DATA, NETCDF4_CLASSIC, and NETCDF4)
#     Default NETCDF4 files use the version 4 disk format (HDF5) and use the new features of the version 4 API.
#     To see how a given file is formatted, you can examine the data_model attribute.
#
#   Modes of opening netcdf files
#     r means read-only; no data can be modified. 
#     w means write; a new file is created, an existing file with the same name is deleted. 
#     a and r+ mean append (in analogy with serial files); an existing file is opened for reading and writing
#
# 
#------------------------------------------------------------------------------------------------------------------
	import netCDF4 as nc
	import numpy as np

	success = False    # switch to True if all steps are completed
	
	fmt = "NETCDF4"
	mode = 'w'
	create_dims = True
	dtype = np.float64
	nvars = len(arrays)
			
	nx = np.size(xvec)
	nz = np.size(zvec)
	nt = 1
	
	# open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(filename, mode, format=fmt)   

	if( create_dims ):
		# Every dimension has a name and length. If we set the dimension length to be 0 or None, 
		# then it takes it as of unlimited size and can grow.
		tdim = nc.createDimension("tdim", nt)    # None triggers unlimited dimension so I can concatenate files across time later
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
		t[:] = tval
		x[:] = xvec[:]
		z[:] = zvec[:]
		
	
	# create the data variable, give it a units attribute and write the data
	for i in range(nvars):
		data = nc.createVariable(varnames[i], dtype, ("tdim","kdim","idim",))
		data.units = dims[i]
		data[:,:] = np.transpose(arrays[i])
	
	
	# print the Dataset object to see what we've got and then close the file
	# print(nc)
	nc.close()
	
	print(filename,' has been created, written to and closed.')
	success = True
	
	return success


def build_compact_matrices_1(N,L):
#-----------------------------------------------------------------------------------------
#  Compute the derivative of f(x) for x in [0,L]  w/o assuming any symmetries at x=0,L
#  Use spectral-like compact scheme of formal 4th order in Lele's paper
#  JOURNAL OF COMPUTATIONAL PHYSICS 103, pp 16-42 (1992) for interior points and
#  for near boundary points use the sixth order schemes from
#  Algorithm 986: A Suite of Compact Finite Difference Schemes
#  ACM Transactions on Mathematical Software, Vol. 44, No. 2, Article 23. Publication date: October 2017
#  MANI MEHRA and KULDIP SINGH PATEL
#-----------------------------------------------------------------------------------------
	import numpy as np
	from scipy.sparse import dia_matrix
	from scipy.linalg import lu_factor
	
	#------------------------------------------------------
	#  the structure of the problem is
	#     A f' = B f
	#  here we build A and B and factor A=LU (w/ pivoting)
	#------------------------------------------------------
	
	#---------------------------------------------------------
	#  parameters for 1st derivative at interior points
	#  (Eqns 2.1 and spectral-like pentadiagonal scheme 3.16)
	#--------------------------------------------------------- 
	alpha = 0.5771439  ; beta = 0.0896406
	a = 1.3025166 ; b = 0.9935500 ; c = 0.0375024 
	
	#------------------------------------------------
	# build A[row,col]  indices start at zero
	#------------------------------------------------
	dx=L/(N-1.); h=dx
	d = np.ones(N)
	data = np.array([beta*d,alpha*d,d,alpha*d,beta*d])
	offsets = np.array([-2,-1,0,1,2])
	A = dia_matrix((data, offsets), shape=(N,N)).toarray()
		
	
	#------------------------------------------------
	# build A[row,col]  indices start at zero
	#------------------------------------------------ 
	data = np.array( [(-c/6.)*d, (-b/4.)*d, (-a/2.)*d, (0.)*d, (a/2.)*d, (b/4.)*d, (c/6.)*d] )/dx
	offsets = np.array([-3,-2,-1,0,1,2,3])
	B = dia_matrix((data, offsets), shape=(N,N)).toarray()
	
	
	#-----------------------------------------------------------------
	# adjust 1st 3 rows (i=0,1,2) of A & B for nonperiodic boundaries
	# 6th order, use results in ACM17 article p. 21/31
	# (fixed typo in matrix summary 1st eqn b=-5/12 not 5/12)
	#-----------------------------------------------------------------
	row=1 ; i=row-1
	alpha = 5. ; a = (-197./60) ; b=(-5./12.) ; c = (5.) ; d = (-5./3.) ; e = (5./12.) ; f = (-1./20)
	A[i,:]=0. ; A[i,i:i+2] = [1.,alpha]
	B[i,:]=0. ; B[i,0:6] = [a/h,b/h,c/h,d/h,e/h,f/h]
	# make the corresponding adjustments for end rows
	A[N-1-i,:] =  np.flip(A[i,:])
	B[N-1-i,:] = -np.flip(B[i,:])     #  note negative sign
	
	row=2 ; i=row-1
	alpha = 2./11. ; a = (-20./33.) ; b=(-35./132.) ; c = (34./33.) ; d = (-7./33) ; e = (2./33.) ; f = (-1./132.)
	A[i,:]=0. ; A[i,i-1:i+2] = [alpha,1.,alpha]
	B[i,:]=0. ; B[i,0:6] = [a/h,b/h,c/h,d/h,e/h,f/h]
	# make the corresponding adjustments for end rows
	A[N-1-i,:] =  np.flip(A[i,:])
	B[N-1-i,:] = -np.flip(B[i,:])     #  note negative sign
	
	row=3 ; i=row-1
	alpha = 1./3. ; a = (-1./36.) ; b=(-14./18.) ; c = (0.) ; d = (14./18.) ; e = (1./36.) 
	A[i,:]=0. ; A[i,i-1:i+2] = [alpha,1.,alpha]
	B[i,:]=0. ; B[i,0:5] = [a/h,b/h,c/h,d/h,e/h]
	# make the corresponding adjustments for end rows
	A[N-1-i,:] =  np.flip(A[i,:])
	B[N-1-i,:] = -np.flip(B[i,:])     #  note negative sign
	
			
	# decompose A into LU
	lu, piv = lu_factor(A)
	
	return lu,piv,B
		
def build_compact_matrices_2(N,L):
#-----------------------------------------------------------------------------------------
#  Compute the derivative of f(x) for x in [0,L]  w/o assuming any symmetries at x=0,L
#  Use spectral-like compact scheme of formal 4th order in Lele's paper
#  JOURNAL OF COMPUTATIONAL PHYSICS 103, pp 16-42 (1992) for interior points and
#  for near boundary points use the sixth order schemes from
#  Algorithm 986: A Suite of Compact Finite Difference Schemes
#  ACM Transactions on Mathematical Software, Vol. 44, No. 2, Article 23. Publication date: October 2017
#  MANI MEHRA and KULDIP SINGH PATEL
#-----------------------------------------------------------------------------------------
	import numpy as np
	from scipy.sparse import dia_matrix
	from scipy.linalg import lu_factor
	
	#------------------------------------------------------
	#  the structure of the problem is
	#     A f' = B f
	#  here we build A and B and factor A=LU (w/ pivoting)
	#------------------------------------------------------
	
	#---------------------------------------------------------
	#  parameters for 2nd derivative at interior points
	#  (Eqns 2.2 and spectral-like pentadiagonal scheme 3.1.12)
	#--------------------------------------------------------- 
	alpha = 0.50209266  ; beta = 0.05569169 
	a =  0.21564935 ; b =  1.7233220 ; c = 0.17659730 
	
	#------------------------------------------------
	# build A[row,col]  indices start at zero
	#------------------------------------------------
	dx=L/(N-1.) ; h2 = dx*dx
	d = np.ones(N)
	data = np.array([beta*d,alpha*d,d,alpha*d,beta*d])
	offsets = np.array([-2,-1,0,1,2])
	A = dia_matrix((data, offsets), shape=(N,N)).toarray()
		
	
	#------------------------------------------------
	# build A[row,col]  indices start at zero
	#------------------------------------------------ 
	data = np.array( [(c/9.)*d, (b/4.)*d, (a)*d, -2.*(c/9. + b/4. + a)*d, (a)*d, (b/4.)*d, (c/9.)*d] )/h2
	offsets = np.array([-3,-2,-1,0,1,2,3])
	B = dia_matrix((data, offsets), shape=(N,N)).toarray()
	
	
	#-----------------------------------------------------------------
	# adjust 1st 3 rows (i=0,1,2) of A & B for nonperiodic boundaries
	# 6th order, use results in ACM17 article p. 24/31
	#-----------------------------------------------------------------
	row=1 ; i=row-1
	alpha = 126./11. ; a = (2077./157.) ; b=(-2943./110.) ; c = (573./44.) ; d = (167./99.) 
	e = (18./11.) ; f = (57./110.) ; g = (-131./1980.)
	A[i,:]=0. ; A[i,i:i+2] = [1.,alpha]
	B[i,:]=0. ; B[i,0:7] = [a/h2,b/h2,c/h2,d/h2,e/h2,f/h2,g/h2]
	# make the corresponding adjustments for end rows
	A[N-1-i,:] =  np.flip(A[i,:])
	B[N-1-i,:] =  np.flip(B[i,:])     #  note absence of negative sign
	
	row=2 ; i=row-1
	alpha = 11./128. ; a = (585./512.) ; b=(-141./64.) ; c = (459./512.) ; d = (9./32.) 
	e = (-81./512.) ; f = (3./64.) ; g = (-3./512)
	A[i,:]=0. ; A[i,i-1:i+2] = [alpha,1.,alpha]
	B[i,:]=0. ; B[i,0:7] = [a/h2,b/h2,c/h2,d/h2,e/h2,f/h2,g/h2]
	# make the corresponding adjustments for end rows
	A[N-1-i,:] =  np.flip(A[i,:])
	B[N-1-i,:] =  np.flip(B[i,:])     #  note absence of negative sign
	
	row=3 ; i=row-1
	alpha = 2./11. ; a = (-3./44.) ; b=(-12./11.) ; c = (-51./22.) ; d = (12./11.) ; e = (3./44.) 
	A[i,:]=0. ; A[i,i-1:i+2] = [alpha,1.,alpha]
	B[i,:]=0. ; B[i,0:5] = [a/h2,b/h2,c/h2,d/h2,e/h2]
	# make the corresponding adjustments for end rows
	A[N-1-i,:] =  np.flip(A[i,:])
	B[N-1-i,:] =  np.flip(B[i,:])     #  note absence of negative sign
	
	# decompose A into LU
	lu, piv = lu_factor(A)
	
	return lu,piv,B	
	
def compact_deriv(f,lu,piv,B):
#-----------------------------------------------------------------------------------------
#  Compute the derivative of f(x) for x in [0,L]  w/o assuming any symmetries at x=0,L
#  Use spectral-like compact scheme of formal 4th order in Lele's paper
#  JOURNAL OF COMPUTATIONAL PHYSICS 103, pp 16-42 (1992)
#    the order of the derivative is determined by the contents of lu and B
#-----------------------------------------------------------------------------------------
	import numpy as np
	from scipy.linalg import lu_solve
	
	# compute the rhs b by matrix multiplication, rhs=Bf
	rhs = B.dot(f)
	
	# solve the linear system for f', store in df
	df = lu_solve((lu, piv), rhs)        # Uses the LU factors.
	    
	return df


def cos_deriv_filtered(f,n,L,frac):
#-----------------------------------------------------------------------------------------
#  compute the nth deriv of even array f assuming equally spaced sampling in [0,L]
#  i.e. f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#  frac = fraction of wavenumber range to filter out  e.g. frac = 0.05
#-----------------------------------------------------------------------------------------
	import numpy as np
	from scipy.fftpack import fft, ifft    
    
	i = 1j                                                   # cmplx(0.,1.)
	N = f.size ;  M = (N-1)*2   # M is array size of even reflection in [0,2L)
	dk = np.pi/L	                                            
	k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
	kmax = np.max(k)  #dk*M/2
	xx = frac*N	
	
	F = np.concatenate([f,  f[1:-1][::-1]])  # F is explicit even extension of data
	FHAT = (i*k)**n * fft(F)                 # nth order deriv operator in Fourier space	    
    
	if( frac > 0. ):
		FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**2))   # taper high k coeffs
	
	df = ifft(FHAT)[0:M/2+1].real    # keep results in [0,L], as real values (imaginary parts = 0)
	return df
    
    
def high_even_deriv(f,order,L,frac):
#-----------------------------------------------------------------------------------------
#  compute the nth deriv of an array f assuming equally spaced sampling in [0,L]
#  f itself does not have zero derivs at 0 and L and so is massaged in a small neighborhood
#  prior to computing the high order (even) derivative. To be used as a dissipation 
#  operator with known inaccuracies near the boundaries. This is acceptable for some
#  applications.
#  frac = fraction of wavenumber range to filter out  e.g. frac = 0.05
#-----------------------------------------------------------------------------------------
	import numpy as np
    
	N = f.size ; F = f.copy()
	nsteps = 3
    
	if(order%2 != 0):
		print('high_even_deriv: requested deriv order must be even. ',n)
		exit() 
    
	if((N-1)%2 != 0):
		print('high_even_deriv: data array must have an odd number of elements on [0,L] ',N,L)
		exit()
    	
	#   smooth the data near end-discontinuity via a diffusion-like process
	kappa = np.zeros_like(f) ; kappa_x = np.zeros_like(f) ; dfdt = np.zeros_like(f)
	h = L/(N-1.); h2=h*h;  gamma=8*h ; kappa_max=1. ; dt = .45*h2/kappa_max
	for i in range(N):
		x = i*h
		if(x <= L/2. ):
			kappa[i] = kappa_max*np.exp(-((x-0.)/gamma)**2) ;  kappa_x[i] = -(2.*(x-0.)/gamma**2)*kappa[i]
		else:
			kappa[i] = kappa_max*np.exp(-((x-L)/gamma)**2) ;  kappa_x[i] = -(2.*(x-L)/gamma**2)*kappa[i]
	

	for n in range(nsteps+1):
		
		i=0
		dfdt[i] = (kappa[i]/h2 - kappa_x[i]/(2*h)) * F[i+1] -(2*kappa[i]/h2)*F[i] + (kappa[i]/h2 + kappa_x[i]/(2*h)) * F[i+1]
		
		i=N-1
		dfdt[i] = (kappa[i]/h2 - kappa_x[i]/(2*h)) * F[i-1] -(2*kappa[i]/h2)*F[i] + (kappa[i]/h2 + kappa_x[i]/(2*h)) * F[i-1]
		
		for i in range(1,N-1):		
			dfdt[i] = (kappa[i]/h2 - kappa_x[i]/(2*h)) * F[i-1] -(2*kappa[i]/h2)*F[i] + (kappa[i]/h2 + kappa_x[i]/(2*h)) * F[i+1]
			
		F[:] = F[:] + dt*dfdt[:]
	
		
	# take high order derivative of this smoothed, cos expandable function
	#df = cos_deriv_filtered(F,order,L,frac)
	df = cos_deriv_filtered(F,order,L,frac)
	
	df = df*(1. - kappa/kappa_max)    # (1 - kappa/kappa_max) is a function that's 1 in the interior and zero at the boundaries
	
	return df,F
