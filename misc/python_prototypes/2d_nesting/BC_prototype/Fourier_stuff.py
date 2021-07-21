#-----------------------------------------------------------------------------------------
#  KBW python scripts for stuff related to Fourier/sin/cos differentiation
#    
#    -----------------------------
#      list of functions defined:
#    -----------------------------
#-----------------------------------------------------------------------------------------

def fourier_wavenumbers(dk,N):
#-----------------------------------------------------------------------------------------
#  construct Fourier wavenumber vector k arranged as [pos 0 negative] wavenumber order
#-----------------------------------------------------------------------------------------
	import numpy as np
	if N%2 == 0:     
		k = dk * np.array([*range(int(N/2+1)),*range(int(-N/2+1),0,1)]) # Python 3
		# k = dk * np.array(range(N/2+1) + range(-N/2+1,0,1))   # Python 2
	else:
		print("fourier_wavenumbers:  N must be even  ",N )
		exit()
	return k

def fourier_filter(k,frac):
#-----------------------------------------------------------------------------------------
#  construct Fourier wavenumber filter that filters out frac*100 % of the
#  largest magnitude Fourier wavenumbers. Assumes [pos 0 negative] wavenumber order.
#  for frac < 0, use the classical 2/3 rule filter
#-----------------------------------------------------------------------------------------
	import numpy as np
	N = k.size
	dk = k[1]-k[0]
	kmax = np.max(k)
	if( frac > 0.0 ):
		xx = frac*N		# fraction of wavenumber range to filter out, e.g. frac=0.05
		filter = (1. - np.exp(-((kmax-np.abs(k))/(xx*dk))**2))	
	elif( frac < 0.0 ):
		filter = np.where(np.abs(k)>(2./3.)*kmax,0.,1)
	elif( frac == 0.0 ):
		filter = np.ones_like(k)
	return filter

	
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
	i = 1j     # cmplx(0.,1.)
    
	#  f(x) can be expanded in a sin series
	if( flag == 1 ):
		F = odd_extend(f)
	#  f(x) can be expanded in a cos series
	if( flag == -1 ):
		F = even_extend(f)
    
	if( flag == 0 ):                   # f(x) can be expanded in a Fourier series                                 
		dk = 2.*np.pi/L
		k = fourier_wavenumbers(dk,N)
		F = f.astype(float)            # was having an endian problem  http://stackoverflow.com/questions/12307429/scipy-fftpack-and-float64
		FHAT = ((i*k)**n) * fft(F)     # have imported scipy's fft,ifft
		df = ifft(FHAT)[0:N].real      # imaginary parts are zero, keep as real array rather than complex with imag=0
		
	if( np.abs(flag) == 1 ):           # f(x) can be expanded in a sine or cosine series
		dk = np.pi/L
		M = (N-1)*2                    # M extended array length
		k = fourier_wavenumbers(dk,M)
		FHAT = ((i*k)**n) * fft(F)     # Fourier transform of nth derivative
		df = ifft(FHAT)[0:M/2+1].real  # nth derivative in [0,L] as real-valued array	
			
	if( np.abs(flag) > 1 ):
		print("dst problem, called with illegal flag value  ",flag)
		exit()
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
	from scipy.fftpack import fft, ifft
    
	N = f.size 
	i = 1j     # cmplx(0.,1.)
    
	#  f(x) can be expanded in a sin series
	if( flag == 1 ):
		F = odd_extend(f)
	#  f(x) can be expanded in a cos series
	if( flag == -1 ):
		F = even_extend(f)
    	
	if(flag == 0):                     # f(x) can be expanded in a Fourier series                                 
		dk = 2.*np.pi/L
		k = fourier_wavenumbers(dk,N)
		F = f.astype(float)            # was having an endian problem  http://stackoverflow.com/questions/12307429/scipy-fftpack-and-float64
		FHAT = ((i*k)**n) * fft(F)     # have imported scipy's fft,ifft
		if(frac>0.): 
			filter = fourier_filter(k,frac)
			FHAT = FHAT*filter
		df = ifft(FHAT)[0:N].real      # imaginary parts are zero, keep as real array rather than complex with imag=0

	if( np.abs(flag) == 1 ):           # f(x) can be expanded in a sine or cosine series
		dk = np.pi/L
		M = (N-1)*2                    # M extended array length
		k = fourier_wavenumbers(dk,M)
		FHAT = ((i*k)**n) * fft(F)     # Fourier transform of nth derivative
		if(frac>0.): 
			filter = fourier_filter(k,frac)
			FHAT = FHAT*filter**n
		#df = ifft(FHAT)[0:M/2+1].real  # nth derivative in [0,L] as real-valued array  #Python 2
		df = ifft(FHAT)[0:int(M/2+1)].real  # nth derivative in [0,L] as real-valued array   #Python 3
		
	if( np.abs(flag) > 1 ):
		print("dst problem, called with illegal flag value  ",flag)
		exit()
	return df


def compute_spectrum(f,flag,L):
#-----------------------------------------------------------------------------------------
#  compute the energy spectrum f_hat * conj(f_hat) (|k|)
#  given f in [0,L) for f periodic or [0,L] for f even or odd symmetry about 0,L
#  i.e. f(0)=f(L) for periodic functions f
#  flag = 0  ==> f is periodic,  f(0)=f(L),  f(L) not explicitly stored
#  flag = 1  ==> f is expandable in sin, f(-x)=f(x) near x=0,L  f(0)=f(L)=0 explicitly stored
#  flag =-1  ==> f is expandable in cos, f(-x)=-f(x) near x=0,L  f(0)=f(L)  explicitly stored
#-----------------------------------------------------------------------------------------
	import numpy as np
	from scipy.fftpack import fft    
    
	if( np.abs(flag) > 1 ):
		print("dst problem, called with illegal flag value  ",flag)
		exit()
	
	N = f.size
	i = 1j                          # cmplx(0.,1.)
    
	if( flag == 1 ):
		F = odd_extend(f)
	if( flag == -1 ):
		F = even_extend(f)
    
	if( flag == 0 ):                #  Fourier series
		dk = 2.*np.pi/L
		M = N                       # M is original array length
		F = f
    
	if( np.abs(flag) == 1 ):        # f(x) can be expanded in a sine or cosine series
		dk = np.pi/L
		M = (N-1)*2                 # M extended array length
		
	k = fourier_wavenumbers(dk,M)     
	FHAT =  fft(F)
			                 
	k = k[0:M/2+1]    								  # keep M/2+1 elements of array, 0, pos k, nyquist
	S = 2 * FHAT[0:M/2+1] * np.conj(FHAT[0:M/2+1])    # multiply by conjugate and double for negative k symmetry
	S[0]=S[0]/2. ; S[-1]=S[-1]/2.                     # zero and Nyquist components don't have negative k counterparts
	S = np.real(S)                                    # store as a real-valued function
	S=S/M                                             # normalize by the array length
	return k,S











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
     print ("2p_diffusion problem, called with illegal flag value")
     print (flag)
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
     print ("compute_spectra problem, called with illegal flag value")
     print (flag)

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
     print ("spectral_filter problem, called with illegal flag value")
     print (flag)
    del F,FHAT,k,dk,M,N
    return df
