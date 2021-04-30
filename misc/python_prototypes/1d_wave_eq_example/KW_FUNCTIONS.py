#-----------------------------------------------------------------------------------------
#  KBW python scripts for general purpose use
#    
#    -----------------------------
#      list of functions defined:
#    -----------------------------
#     dst(f,n,L,flag)                compute deriv of f using spectral transforms
#
#     int_a2b(x,y,a,b,M)             given x,y; compute integral yi(x) dxi from a to b, 
#                                     where yi,xi are interpolated in [a,b] w/ M pts
#-----------------------------------------------------------------------------------------




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
    from scipy.fftpack import fft, ifft, rfft, dct
    
    N = f.size
    i = 1j                                                   # cmplx(0.,1.)
    if flag == 0:                                            #  Fourier series
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
