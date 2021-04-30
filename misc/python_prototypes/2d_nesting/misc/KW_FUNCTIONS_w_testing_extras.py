#-----------------------------------------------------------------------------------------
#  KBW python scripts for general purpose use
#    
#    -----------------------------
#      list of functions defined:
#    -----------------------------
#     df=dst(f,n,L,flag)                                        compute deriv of f using spectral transforms
#     df=dst_filtered(f,n,L,flag,frac)                          dst, but filter a fraction of high wavenumber components
#     f_x,f_z = grad(f,Lx,Lz,frac)                              compute gradient of f[x,z]   (f cos expandable in [0,Lx],[0,Lz] )
#     div = divergence(u,w,Lx,Lz,frac)                          compute div = u_x + w_z    (u,w cos expandable in [0,Lx],[0,Lz] )
#     grad2_f = grad2(f,Lx,Lz,frac)                             compute f_xx + f_zz          (f cos expandable in [0,Lx],[0,Lz] )
#     inversion_matrix = create_inversion_matrix(nx,nz,Lx,Lz)   fill 2d array -1/(k^2 + m^2) for even-extended, periodic data
#     p=poisson_invert_2d(f,L)                                  solve p_xx + p_zz = f  assuming f can be even extended in x,z + homog. Neumann BCs
#
#     other misc test routines
#-----------------------------------------------------------------------------------------

def compute_deriv(f,L,frac):
#----------------------------------------------------------------------------
#    compute f' of f(x) defined on the discrete closed interval [0,L] 
#      N.B.   f is neither purely even nor odd
#             use dst_filtered(f[:,k],order,L,flag,frac)
#----------------------------------------------------------------------------
	import numpy as np
	N = f.size
	NMID = (N-1)/2 + 1
	x = np.linspace(0,L,N) ; dx = x[1]-x[0] ; gamma=dx*(N-1)/4.
	df = np.zeros_like(f) 
	u  = np.zeros_like(f); u_x = np.zeros_like(f) ; u_xx = np.zeros_like(f)
		
	
	#-------------------------------------------------------------------------
	#  (1)  for g=(1-e^u), define u(x) and evaluate its first two derivatives
	#-------------------------------------------------------------------------
	p = 8 ; pm1 = p-1 ; pm2 = p - 2
	
	u[0:NMID] = -( (x[0:NMID]- 0.)/gamma)**p
	u[NMID:N] = -( (x[NMID:N]- L )/gamma)**p
	
	u_x[0:NMID] = -(p/gamma)*( (x[0:NMID]- 0.)/gamma )**pm1
	u_x[NMID:N] = -(p/gamma)*( (x[NMID:N]- L )/gamma )**pm1
	
	u_xx[0:NMID] = -(p/gamma)*(pm1/gamma)*( (x[0:NMID]- 0.)/gamma )**pm2
	u_xx[NMID:N] = -(p/gamma)*(pm1/gamma)*( (x[NMID:N]- L )/gamma )**pm2
	
	#---------------------------------------------------------------
	#  (2)  evaluate g(x) and its first two derivatives
	#       N.B. g=g'=0 at x=0,L
	#---------------------------------------------------------------
	g = 1.0 - np.exp(u)
	g_x = -u_x * np.exp(u)
	g_xx = -( u_x**2 + u_xx ) * np.exp(u)
	
	#----------------------------------------------------------------------
	#  (3)  form the modified function F(x) = g(x)*f(x)
	#       F = 0 at x=0,L
	#       F' = 0 at x=0,L  ==> F can be cos expanded and differentiated
	#       F" = fg" + 2f'g' + gf" = 0 at x=0,L  when g"=0 at ends, i.e. p=4
	#       F"' = fg"' + gf"' + 3g"f' + 3 f'g" = fg"' at x=0,L when p=4, continuous
	#
	#       in the interior,
	#       F' = gf' + fg'  ==>  f' = (F'- fg')/g  except where g=0
	#----------------------------------------------------------------------	
	F = g[:]*f[:]
	
	
	#---------------------------------------------------------------
	#  (4) take FFT derivative of even extension of F
	#      filter w/ "frac", keep results in [0,L]
	#---------------------------------------------------------------	
	flag = -1  #  use FFT of even expansion of F
	order = 1  #  compute first derivative
	F_x  = dst_filtered(F,order,L,flag,frac)
		
	
	#---------------------------------------------------------------
	#  (5)  extract f' at all interior gridpoints
	#---------------------------------------------------------------
	#df[1:N-2] = ( F_x[1:N-2] - f[1:N-2]*g_x[1:N-2] )/g[1:N-2]
	df[2:N-3] = ( F_x[2:N-3] - f[2:N-3]*g_x[2:N-3] )/g[2:N-3]
	
	
	#---------------------------------------------------------------
	#  (6)  use one-sided finite difference estimates for end values
	#       replace a few near-boundary values w/ FD estimates
	#       based on interior information only
	#---------------------------------------------------------------
	# end points
	df[0] = (-3.*f[0] + 4.*f[1] - f[2])/(2.*dx)
	df[N-1] = (3.*f[N-1] - 4.*f[N-2] + f[N-3])/(2.*dx)
	
	# 1st interior points
	i=1
	df[i] = (f[i+1]-f[i-1])/(2.*dx)
	i=N-2
	df[i] = (f[i+1]-f[i-1])/(2.*dx)
	
	npts = np.int(.5*gamma/dx)
	for i in range(2,2+npts):
		df[i] = (f[i-2] - 8.*f[i-1] + 8.*f[i+1] - f[i+2])/(12.*dx)
		
	for i in range(N-2-npts,N-2):
		df[i] = (f[i-2] - 8.*f[i-1] + 8.*f[i+1] - f[i+2])/(12.*dx)	
			
	return df     
	
	
	
	
def grad(f,Lx,Lz,frac):
#-------------------------------------------------------------------------------------------------
#    compute grad f[x,z] = f_x and f_z using filtered dst routines w/ filtering parameter frac
#    f[x,z] w/ no particular symmetry on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(f)
	flag = -1       #  use even expansions for FFTs
	order = 1       #  compute first derivatives
	f_x = np.zeros_like(f) ; f_z = np.zeros_like(f)
	
	# compute f_x
	for k in range(nz):
		#f_x[:,k] = dst_filtered(f[:,k],order,Lx,flag,frac)
		f_x[:,k] = compute_deriv(f[:,k],Lx,frac)
	
	# compute f_z
	for i in range(nx):
		#f_z[i,:] = dst_filtered(f[i,:],order,Lz,flag,frac)
		f_z[i,:] = compute_deriv(f[i,:],Lz,frac)			
	return f_x,f_z

	
def divergence(u,w,Lx,Lz,frac):
#-------------------------------------------------------------------------------------------------
#    compute div = u_x + w_z using filtered dst routines w/ filtering parameter frac
#    u[x,z] and v[x,z] w/ no particular symmetry on [0,Lx] [0,Lz]
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(u)
	flag = -1       #  use even expansions for FFTs
	order = 1       #  compute first derivatives
	u_x = np.zeros_like(u) ; w_z = np.zeros_like(w)
	
	# compute u_x
	for k in range(nz):
		#u_x[:,k] = dst_filtered(u[:,k],order,Lx,flag,frac)
		u_x[:,k] = compute_deriv(u[:,k],Lx,frac)
	
	# compute f_z
	for i in range(nx):
		#w_z[i,:] = dst_filtered(w[i,:],order,Lz,flag,frac)
		w_z[i,:] = compute_deriv(w[i,:],Lz,frac)
	
	div = u_x + w_z	
	del u_x,w_z	
	return div


def grad2(f,Lx,Lz,frac):
#-------------------------------------------------------------------------------------------------
#    compute grad^2 f[x,z] = f_xx + f_zz using filtered dst routines w/ filtering parameter frac
#-------------------------------------------------------------------------------------------------
	import numpy as np
	nx,nz = np.shape(f)
	fx,fz = grad(f,Lx,Lz,frac)
	grad2_f = divergence(fx,fz,Lx,Lz,frac)
	return grad2_f


def compute_derivs(f,L,frac):
#----------------------------------------------------------------------------
#    compute f' and f" of f(x) defined on the discrete closed interval [0,L] 
#      N.B.   f is neither purely even nor odd
#             use dst_filtered(f[:,k],order,L,flag,frac)
#----------------------------------------------------------------------------
	import numpy as np
	N = f.size
	NMID = (N-1)/2 + 1
	x = np.linspace(0,L,N) ; dx = x[1]-x[0] ; gamma=dx*(N-1)/4.
	df = np.zeros_like(f) ; ddf = np.zeros_like(f)
	u  = np.zeros_like(f); u_x = np.zeros_like(f) ; u_xx = np.zeros_like(f)
		
	
	#-------------------------------------------------------------------------
	#  (1)  for g=(1-e^u), define u(x) and evaluate its first two derivatives
	#-------------------------------------------------------------------------
	p = 8 ; pm1 = p-1 ; pm2 = p - 2
	
	u[0:NMID] = -( (x[0:NMID]- 0.)/gamma)**p
	u[NMID:N] = -( (x[NMID:N]- L )/gamma)**p
	
	u_x[0:NMID] = -(p/gamma)*( (x[0:NMID]- 0.)/gamma )**pm1
	u_x[NMID:N] = -(p/gamma)*( (x[NMID:N]- L )/gamma )**pm1
	
	u_xx[0:NMID] = -(p/gamma)*(pm1/gamma)*( (x[0:NMID]- 0.)/gamma )**pm2
	u_xx[NMID:N] = -(p/gamma)*(pm1/gamma)*( (x[NMID:N]- L )/gamma )**pm2
	
	#---------------------------------------------------------------
	#  (2)  evaluate g(x) and its first two derivatives
	#       N.B. g=g'=0 at x=0,L
	#---------------------------------------------------------------
	g = 1.0 - np.exp(u)
	g_x = -u_x * np.exp(u)
	g_xx = -( u_x**2 + u_xx ) * np.exp(u)
	
	#----------------------------------------------------------------------
	#  (3)  form the modified function F(x) = g(x)*f(x)
	#       F = 0 at x=0,L
	#       F' = 0 at x=0,L  ==> F can be cos expanded and differentiated
	#       F" = fg" + 2f'g' + gf" = 0 at x=0,L  when g"=0 at ends, i.e. p=4
	#       F"' = fg"' + gf"' + 3g"f' + 3 f'g" = fg"' at x=0,L when p=4, continuous
	#
	#       in the interior,
	#       F' = gf' + fg'  ==>  f' = (F'- fg')/g  except where g=0
	#----------------------------------------------------------------------	
	F = g[:]*f[:]
	
	
	#---------------------------------------------------------------
	#  (4) take FFT derivative of even extension of F
	#      filter w/ "frac", keep results in [0,L]
	#---------------------------------------------------------------	
	flag = -1  #  use FFT of even expansion of F
	order = 1  #  compute first derivative
	F_x  = dst_filtered(F,order,L,flag,frac)
	
	order = 2  #  compute second derivative
	F_xx = dst_filtered(F,order,L,flag,frac)
	
	
	#---------------------------------------------------------------
	#  (5)  extract f' at all interior gridpoints
	#---------------------------------------------------------------
	#df[1:N-2] = ( F_x[1:N-2] - f[1:N-2]*g_x[1:N-2] )/g[1:N-2]
	df[2:N-3] = ( F_x[2:N-3] - f[2:N-3]*g_x[2:N-3] )/g[2:N-3]
	
	
	#---------------------------------------------------------------
	#  (6)  use one-sided finite difference estimates for end values
	#       replace a few near-boundary values w/ FD estimates
	#       based on interior information only
	#---------------------------------------------------------------
	# end points
	df[0] = (-3.*f[0] + 4.*f[1] - f[2])/(2.*dx)
	df[N-1] = (3.*f[N-1] - 4.*f[N-2] + f[N-3])/(2.*dx)
	
	# 1st interior points
	i=1
	df[i] = (f[i+1]-f[i-1])/(2.*dx)
	i=N-2
	df[i] = (f[i+1]-f[i-1])/(2.*dx)
	
	npts = np.int(.5*gamma/dx)
	for i in range(2,2+npts):
		df[i] = (f[i-2] - 8.*f[i-1] + 8.*f[i+1] - f[i+2])/(12.*dx)
		
	for i in range(N-2-npts,N-2):
		df[i] = (f[i-2] - 8.*f[i-1] + 8.*f[i+1] - f[i+2])/(12.*dx)
	
	
	
	#---------------------------------------------------------------
	#  (7)  extract f'' at all interior gridpoints
	#---------------------------------------------------------------
	#df[1:N-2] = ( F_xx[1:N-2] - f[1:N-2]*g_xx[1:N-2] - 2.*df[1:N-2]*g_x[1:N-2] )/g[1:N-2]
	ddf[2:N-3] = ( F_xx[2:N-3] - f[2:N-3]*g_xx[2:N-3] - 2.*df[2:N-3]*g_x[2:N-3])/g[2:N-3]
	
	# end points
	ddf[0] = (2*f[0] -5.*f[1] + 4.*f[2] -1.*f[3])/(dx**2)
	ddf[N-1] = (2*f[N-1] -5.*f[N-2] + 4.*f[N-3] -1.*f[N-4])/(dx**2)
	
	
	# 1st interior points
	i=1
	ddf[i] = (f[i+1] -2.*f[i] + f[i-1])/(dx**2)
	i=N-2
	ddf[i] = (f[i+1] -2.*f[i] + f[i-1])/(dx**2)
	
	npts = np.int(.5*gamma/dx)
	for i in range(2,2+npts):
		ddf[i] = (-1.*f[i-2] +16.*f[i-1] - 30.*f[i] +16.*f[i+1] -1.*f[i+2])/(12.*dx**2)
		
	for i in range(N-2-npts,N-2):
		ddf[i] = (-1.*f[i-2] +16.*f[i-1] - 30.*f[i] +16.*f[i+1] -1.*f[i+2])/(12.*dx**2)
	
			
	return df,ddf  


def dst(f,n,L,flag):
#-----------------------------------------------------------------------------------------
#  compute the nth deriv of array f assuming equally spaced
#  sampling with periodicity/even/odd distance equal to L
#  i.e. f(0)=f(L) for periodic functions f
#  flag = 0  ==> f is periodic,  f(0)=f(L),  f(L) not explicitly stored
#  flag = 1  ==> f is expandable in sin, f(-x)=f(x) near x=0,L  f(0)=f(L)=0 explicitly stored
#  flag =-1  ==> f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft
    
    N = f.size
    i = 1j                                                   # cmplx(0.,1.)
    if flag == 0:                                            # Fourier series
     dk = 2.*np.pi/L
     if N%2 == 0:
      k = dk * np.array(range(N/2+1) + range(-N/2+1,0,1))    #  N even, usual case
     else:
      k = dk * np.array(range(1,(N-1)/2+1,1) + range(-(N-1)/2,1,1))  # N odd (??)
     F = f.astype(float)                       #  was having an endian problem  http://stackoverflow.com/questions/12307429/scipy-fftpack-and-float64
     FHAT = ((i*k)**n) * fft(F)                # have imported scipy's fft,ifft
     df = ifft(FHAT)[0:N].real
    elif flag == 1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
     F = np.concatenate([f, -f[1:-1][::-1]])                 # F is odd extension of data
     FHAT = ((i*k)**n) * fft(F)                              # have imported scipy's fft,ifft
     df = ifft(FHAT)[0:M/2+1].real 
    elif flag == -1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
     F = np.concatenate([f,  f[1:-1][::-1]])                 # F is even extension of data
     FHAT = ((i*k)**n) * fft(F)                              # have imported scipy's fft,ifft
     df = ifft(FHAT)[0:M/2+1].real
    else:
     print "dst problem, called with illegal flag value"
     print flag
    del F,FHAT,k
    return df

def dst_filtered(f,n,L,flag,frac):
#-----------------------------------------------------------------------------------------
#  compute the nth deriv of array f assuming equally spaced
#  sampling with periodicity/even/odd distance equal to L
#  i.e. f(0)=f(L) for periodic functions f
#  flag = 0  ==> f is periodic,  f(0)=f(L),  f(L) not explicitly stored
#  flag = 1  ==> f is expandable in sin, f(-x)=f(x) near x=0,L  f(0)=f(L)=0 explicitly stored
#  flag =-1  ==> f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft, rfft, dct    
    
    N = f.size
    xx = frac*N		# fraction of wavenumber range to filter out 0.05
    
    i = 1j                                                   # cmplx(0.,1.)
    
    if flag == 0:                                            #  Fourier series
     dk = 2.*np.pi/L
     if N%2 == 0:
      k = dk * np.array(range(N/2+1) + range(-N/2+1,0,1))    #  N even, usual case
     else:
      k = dk * np.array(range(1,(N-1)/2+1,1) + range(-(N-1)/2,1,1))  # N odd (??) Fourier should always be even
     
     kmax = np.max(k)
     FHAT = ((i*k)**n) * fft(f)                              # have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**2))
     df = ifft(FHAT)[0:N].real
    elif flag == 1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
     
     kmax = np.max(k)
     F = np.concatenate([f, -f[1:-1][::-1]])                 # F is odd extension of data
     FHAT = ((i*k)**n) * fft(F)                              # have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**2))
     df = ifft(FHAT)[0:M/2+1].real 
    elif flag == -1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
    
     kmax = np.max(k)
     F = np.concatenate([f,  f[1:-1][::-1]])                 # F is even extension of data
     FHAT = ((i*k)**n) * fft(F)                              # have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**2))
     df = ifft(FHAT)[0:M/2+1].real
    else:
     print "dst_filtered problem, called with illegal flag value"
     print flag
    del F,FHAT,k,dk,M,N
    return df


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
	print(nc)
	nc.close()
	
	print(filename,' has been created, written to and closed.')
	success = True
	
	return success


def poisson_invert(f,L):
#-----------------------------------------------------------------------------------------
#  solve p_xx = f(x)   p_x(0)=p_x(L)=0, equally spaced closed interval grid p,f, even
#  flag =-1  ==> f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft, rfft, dct
    
    N = f.size
    i = 1j                                                   # cmplx(0.,1.)
    dk = np.pi/L
    M = (N-1)*2                                              # extended array length
    k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
    F = np.concatenate([f,  f[1:-1][::-1]])                  # F is even extension of data
    FHAT = fft(F)                                            # FT of even extended periodic data vector
    phat = np.zeros_like(FHAT)
    
    for i in range(1,M):       # skip 0 wavenumber location
    	phat[i] = -(1./k[i]**2) * FHAT[i]                    # integrate twice wrt x in wavenumber space
    
    
	p = ifft(phat)[0:M/2+1].real                             # inverse transform, keep cos series
	p_ext = ifft(phat).real
    return p,p_ext,F,FHAT,k



def build_compact_matrices_1(N,L):
#-----------------------------------------------------------------------------------------
#  Compute the derivative of f(x) for x in [0,L]  w/o assuming any symmetries at x=0,L
#  Use spectral-like compact scheme of formal 4th order in Lele's paper
#  JOURNAL OF COMPUTATIONAL PHYSICS 103, pp 16-42 (1992)
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
	# 6th order one sided, eqns 4.1.1 and 4.1.4
	#-----------------------------------------------------------------
	alpha = 3.0 ; a = (-17./6.) ; b=(3./2.) ; c=b ; d = (-1./6.) 
	for row in range(1,4):
		i=row-1   	         # row = eqn number, i 0-based index	 
		A[i,:]=0. ; A[i,i:i+2] = [1,alpha]
		B[i,:]=0. ; B[i,i:i+4] = [a/h,b/h,c/h,d/h]
		# make the corresponding adjustments to the last bottom rows
		A[N-1-i,:] =  np.flip(A[i,:])
		B[N-1-i,:] = -np.flip(B[i,:])     #  note negative sign
	
	
	# decompose A into LU
	lu, piv = lu_factor(A)
	
	return lu,piv,B
		
def build_compact_matrices_2(N,L):
#-----------------------------------------------------------------------------------------
#  Compute the 2nd derivative of f(x) for x in [0,L]  w/o assuming any symmetries at x=0,L
#  Use spectral-like compact scheme of formal 4th order in Lele's paper
#  JOURNAL OF COMPUTATIONAL PHYSICS 103, pp 16-42 (1992)
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
	# 3rd order one sided, eqns 4.3.1 
	#-----------------------------------------------------------------
	alpha = 11. ; a = 13. ; b = -27. ; c = 15.; d = -1.
	for row in range(1,4):
		i=row-1   	         # row = eqn number, i 0-based index	 
		A[i,:]=0. ; A[i,i:i+2] = [1,alpha]
		B[i,:]=0. ; B[i,i:i+4] = [a/h2,b/h2,c/h2,d/h2]
		# make the corresponding adjustments to the last bottom rows
		A[N-1-i,:] =  np.flip(A[i,:])
		B[N-1-i,:] =  np.flip(B[i,:])     #  no negative sign
	
	
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



def ddx_terrible(f,L,frac):
#-----------------------------------------------------------------------------------------
#  Compute the derivative of f(x) for x in [0,L]  w/o assuming any symmetries at x=0,L
#  Use the formal definition of the cos expansion coefficients for f' and integrate
#  by parts twice. This yields some boundary terms involving f and the expansion
#  coeffs of F(x), the indefinite integral of f.
#-----------------------------------------------------------------------------------------
	import numpy as np
	from scipy.fftpack import fft, ifft 
	from scipy.integrate import cumtrapz  
    
	i = 1j          # cmplx(0.,1.)
	pi = np.pi
	N = f.size ; x=np.linspace(0,L,N)
	xx = frac*N		# fraction of wavenumber range to filter out 0.05
    
	dk = pi/L
	M = (N-1)*2                                             # extended array length
	k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
	kmax = np.max(k)
	
	# simple integration scheme
	F = cumtrapz(f,x,initial=0.0)  # Cumulatively integrate f(x) using the composite trapezoidal rule.
	
	# form EVEN extension of F 
	F = np.concatenate([F, F[1:-1][::-1]])    # F is even extension of integrated data
	
	# get k^2*cos series expansion coeffs via (ik)**2 * Fourier coeffs
	FHAT = (i*k)**2 * fft(F)
	
	# add in the boundary terms from integration by parts
	alpha = 2./L
	for ii in range(M):
		FHAT[ii] = alpha*( f[-1]*np.cos(ii*pi) - f[0] ) + FHAT[ii]
     
     
    #  filter if desired
	if( frac > 0. ):
		FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**2))
    
    # inverse FFT, keep results in [0,L] as real values (imaginary parts are zero)
	df = ifft( FHAT )[0:M/2+1].real 
    
	return df






def smooth_near_bdry(f,nwidth):
#-----------------------------------------------------------------------------------------
#  apply near boundary smoothing to f, smooth over ! "nwidth" grid points
#-----------------------------------------------------------------------------------------
	import numpy as np
	n = f.size    # assume odd because we are using cos/sin expansions
	nb2 = n/2
	npts = 4*nwidth
	
	f_smooth = f.copy()
	for i in range(1,npts):
		alpha = np.exp(-((i-0)/nwidth)**2)     #  1 at bdry,  --> 0 in far interior
		beta = 1.0 - alpha                     #  0 at bdry,  --> 1 in far interior
		f_smooth[i] = (alpha*f[i-1] + beta*f[i] + alpha*f[i+1])/(2.*alpha+beta)
	
	for i in range(n-npts,n-1):
		alpha = np.exp(-((i-(n-1))/nwidth)**2)     #  1 at bdry,  --> 0 in far interior
		beta = 1.0 - alpha                         #  0 at bdry,  --> 1 in far interior
		f_smooth[i] = (alpha*f[i-1] + beta*f[i] + alpha*f[i+1])/(2.*alpha+beta)
	
	return f_smooth






def integrate_a2b(x,y,a,b,M):
    import numpy as np
    from scipy import interpolate
    from scipy.integrate import simps,trapz
    f = interpolate.interp1d(x,y,kind='cubic')    # create an interpolator
    xi = np.linspace(a,b,M)                       # size M grid that includes endpts
    yi = f(xi)                                    # interpolated values
    I = simps(yi,xi)
    del f,xi,yi
    return I 


def diffusion_2p(f,p,L,flag,nu_star,frac):
#-----------------------------------------------------------------------------------------
#  compute sign * nu* times the 2pth deriv of array f assuming equally spaced
#  sampling with periodicity/even/odd distance equal to L
#  i.e. f(0)=f(L) for periodic functions f
#  flag = 0  ==> f is periodic,  f(0)=f(L),  f(L) not explicitly stored
#  flag = 1  ==> f is expandable in sin, f(-x)=f(x) near x=0,L  f(0)=f(L)=0 explicitly stored
#  flag =-1  ==> f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft, rfft, dct
        
    n = 2*p													 # even order of derivative operator
    sign = (-1)**(p-1)										 # sign of diffusive term p=1==>diffusion ==> +
    														 #						  p=2==>biharmonic==> -  etc
    														 
    zz = nu_star**(1.0/n)			# i.e. xx**n = nu_star   absorb this into "k**n" to avoid overflow prior to multiplying by small nu_star
    								#  i.e. k**n  will really be nu_star * k**n in what follows below 
    
    N = f.size
    xx = frac*N		# fraction of wavenumber range to filter out 0.05
    
    i = 1j                                                   # cmplx(0.,1.)
    
    if flag == 0:                                            #  Fourier series
     dk = zz*2.*np.pi/L				#  NB zz factor
     if N%2 == 0:
      k = dk * np.array(range(N/2+1) + range(-N/2+1,0,1))    #  N even, usual case
     else:
      k = dk * np.array(range(1,(N-1)/2+1,1) + range(-(N-1)/2,1,1))  # N odd (??)
     kmax = np.max(k)
     F = f.astype(float)                    #  was having an endian problem  http://stackoverflow.com/questions/12307429/scipy-fftpack-and-float64
     FHAT = fft(F)							# have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**(2*p)))   # 2p filter
     FHAT = ((i*k)**n) * FHAT                #  mult by i v_star k^2p
     diff = sign * ifft(FHAT)[0:N].real
    elif flag == 1:
     dk = zz*np.pi/L				#  NB zz factor
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
     kmax = np.max(k)
     F = np.concatenate([f, -f[1:-1][::-1]])                 # F is odd extension of data
     FHAT = fft(F)							# have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**(2*p)))   # 2p filter
     FHAT = ((i*k)**n) * FHAT                #  mult by i v_star k^2p
     diff = sign * ifft(FHAT)[0:M/2+1].real 
    elif flag == -1:
     dk = zz*np.pi/L				#  NB zz factor
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
     kmax = np.max(k)
     F = np.concatenate([f,  f[1:-1][::-1]])                 # F is even extension of data
     FHAT = fft(F)							# have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**(2*p)))   # 2p filter
     FHAT = ((i*k)**n) * FHAT                #  mult by i v_star k^2p
     diff = sign * ifft(FHAT)[0:M/2+1].real
    else:
     print "2p_diffusion problem, called with illegal flag value"
     print flag
    del F,FHAT,k
    return diff


def compute_spectra(f,flag,L,nu_star,p):
#-----------------------------------------------------------------------------------------
#  compute the nth deriv of array f assuming equally spaced
#  sampling with periodicity/even/odd distance equal to L
#  i.e. f(0)=f(L) for periodic functions f
#  flag = 0  ==> f is periodic,  f(0)=f(L),  f(L) not explicitly stored
#  flag = 1  ==> f is expandable in sin, f(-x)=f(x) near x=0,L  f(0)=f(L)=0 explicitly stored
#  flag =-1  ==> f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft, rfft, dct    
    
    N = f.size
    i = 1j                                                   # cmplx(0.,1.)
    
    if flag == 0:                                            #  Fourier series
     dk = 2.*np.pi/L
     if N%2 == 0:
      k = dk * np.array(range(N/2+1) + range(-N/2+1,0,1))    #  N even, usual case
     else:
      k = dk * np.array(range(1,(N-1)/2+1,1) + range(-(N-1)/2,1,1))  # N odd (??) Fourier should always be even     
     FHAT =  fft(f)                              			# have imported scipy's fft,ifft
    elif flag == 1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )     
     F = np.concatenate([f, -f[1:-1][::-1]])                 # F is odd extension of data
     FHAT =  fft(F)                              			# have imported scipy's fft,ifft
    elif flag == -1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
     F = np.concatenate([f,  f[1:-1][::-1]])                 # F is even extension of data
     FHAT =  fft(F)                              			# have imported scipy's fft,ifft
    else:
     print "compute_spectra problem, called with illegal flag value"
     print flag

    k = k[0:M/2+1]    								# keep M/2+1 elements of array, 0, pos k, nyquist
    E = 2 * FHAT[0:M/2+1] * np.conj(FHAT[0:M/2+1])
    E[0]=E[0]/2. ; E[-1]=E[-1]/2. ; E = np.real(E) ; E=E/N
    D = nu_star * k[0:M/2+2]**(2*p)
    return k,E,D


def spectral_filter(f,frac,L,flag,p):
#-----------------------------------------------------------------------------------------
#  compute the nth deriv of array f assuming equally spaced
#  sampling with periodicity/even/odd distance equal to L
#  i.e. f(0)=f(L) for periodic functions f
#  flag = 0  ==> f is periodic,  f(0)=f(L),  f(L) not explicitly stored
#  flag = 1  ==> f is expandable in sin, f(-x)=f(x) near x=0,L  f(0)=f(L)=0 explicitly stored
#  flag =-1  ==> f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft, rfft, dct    
    
    N = f.size
    xx = frac*N    #  is fraction of wavenumber range to filter out, e.g. 0.05
    
    i = 1j                                                   # cmplx(0.,1.)
    
    if flag == 0:                                            #  Fourier series
     dk = 2.*np.pi/L
     if N%2 == 0:
      k = dk * np.array(range(N/2+1) + range(-N/2+1,0,1))    #  N even, usual case
     else:
      k = dk * np.array(range(1,(N-1)/2+1,1) + range(-(N-1)/2,1,1))  # N odd (??) Fourier should always be even
     
     kmax = np.max(k)
     FHAT = fft(f)                              # have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**(2*p)))
     df = ifft(FHAT)[0:N].real
    elif flag == 1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
     
     kmax = np.max(k)
     F = np.concatenate([f, -f[1:-1][::-1]])                 # F is odd extension of data
     FHAT = fft(F)                              # have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**(2*p)))
     df = ifft(FHAT)[0:M/2+1].real 
    elif flag == -1:
     dk = np.pi/L
     M = (N-1)*2                                             # extended array length
     k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
    
     kmax = np.max(k)
     F = np.concatenate([f,  f[1:-1][::-1]])                 # F is even extension of data
     FHAT = fft(F)                              # have imported scipy's fft,ifft
     FHAT = FHAT * (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**(2*p)))
     df = ifft(FHAT)[0:M/2+1].real
    else:
     print "spectral_filter problem, called with illegal flag value"
     print flag
    del F,FHAT,k,dk,M,N
    return df

    
    
    
