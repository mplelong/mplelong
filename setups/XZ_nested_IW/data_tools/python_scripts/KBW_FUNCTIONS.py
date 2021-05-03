#-----------------------------------------------------------------------------------------
#  KBW python scripts for general purpose use
#    this file needs to be in the path defined by env variable PYTHONPATH
#    
#    -----------------------------
#      list of functions defined:
#    -----------------------------
#     vort_z(u,v,Lx,Ly):             compute dvdx - dudy                         
#
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
     print ("dst problem, called with illegal flag value")
     print (flag)
    del F,FHAT,k
    return df

#---------------------------------------------------------------------------------------
#  compute z vorticity
#---------------------------------------------------------------------------------------
def vort_z(u,v,Lx,Ly):     #  dvdx - dudy
#  compute vertical vorticity using spectral differentiation given planar data
    import os,sys
    import numpy as np        
    [ny,nx]=v.shape     # input array shape      
#-----------------------------------------------------------------------------------------
#   set the derivative flags for this problem:
#   for even data flag=-1, for odd data flag=1, for periodic data flag=0
#-----------------------------------------------------------------------------------------
    flagvx=0;             #  all fields periodic in x
    flaguy=-1;            #  u even/cos in y
#-----------------------------------------------------------------------------------------
#   do the y derivatives first
#-----------------------------------------------------------------------------------------
    order=1                      # first derivatives
    duy=0*u.copy()
    for i in range(nx):
        duy[:,i] = dst( u[:,i],order,Ly,flaguy )   
#-----------------------------------------------------------------------------------------
#   do the x derivatives next
#-----------------------------------------------------------------------------------------
    dvx=0*v.copy(); 
    order=1                      # first derivatives
    for j in range(ny):
        dvx[j,:] = dst( v[j,:],order,Lx,flagvx )    

    ans = dvx - duy         
    return ans


def vort_z_fd(u,v,Lx,Ly):     #  dvdx - dudy
#  compute vertical vorticity using spectral differentiation given planar data
    import os,sys
    import numpy as np        
    [ny,nx]=v.shape     # input array shape      
#-----------------------------------------------------------------------------------------
#   compute grid spacing (assume closed intervals, cos expansions)
#-----------------------------------------------------------------------------------------
    dx=Lx/(nx-1.)
    dy=Ly/(ny-1.)
    
#-----------------------------------------------------------------------------------------
#   do the y derivatives first
#-----------------------------------------------------------------------------------------
    order=1                      # first derivatives
    duy=np.zeros_like(u) 
    for i in range(nx):
        duy[:,i] = deriv_fd( u[:,i],dy )   
#-----------------------------------------------------------------------------------------
#   do the x derivatives next
#-----------------------------------------------------------------------------------------
    dvx=np.zeros_like(v) 
    order=1                      # first derivatives
    for j in range(ny):
        dvx[j,:] = deriv_fd( v[j,:],dx )    

    ans = dvx - duy         
    return ans


#---------------------------------------------------------------------------------------
#  compute 2d u dot grad phi  u*phi_x + v*phi_y
#---------------------------------------------------------------------------------------
def u_dot_grad_2d(u,v,phi,Lx,Ly,xflag,yflag):    
#  compute horizontal u dot grad phi
    import os,sys
    import numpy as np 
           
    [ny,nx]=v.shape                    # input array shape
    deriv_order=1                      # first derivatives 
         
#-----------------------------------------------------------------------------------------
#   for even data flag=-1, for odd data flag=1, for periodic data flag=0
#-----------------------------------------------------------------------------------------
    if(yflag=='c'): flag=-1
    if(yflag=='s'): flag= 1
    if(yflag=='p'): flag= 0
	
#-----------------------------------------------------------------------------------------
#   do the y derivatives first
#-----------------------------------------------------------------------------------------
    phi_y=np.zeros_like(phi)
    for i in range(nx):
        phi_y[:,i] = dst( phi[:,i],deriv_order,Ly,flag ) 
        
#-----------------------------------------------------------------------------------------
#   for even data flag=-1, for odd data flag=1, for periodic data flag=0
#-----------------------------------------------------------------------------------------
    if(xflag=='c'): flag=-1
    if(xflag=='s'): flag= 1
    if(xflag=='p'): flag= 0
          
#-----------------------------------------------------------------------------------------
#   do the x derivatives next
#-----------------------------------------------------------------------------------------
    phi_x=np.zeros_like(phi) 
    for j in range(ny):
        phi_x[j,:] = dst( phi[j,:],deriv_order,Lx,flag )    

#-----------------------------------------------------------------------------------------
#   multiply and add results, clean up
#-----------------------------------------------------------------------------------------
    ans = u*phi_x + v*phi_y
    del phi_x,phi_y,flag,deriv_order         
    return ans









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

def deriv_fd(f,h):
#-----------------------------------------------------------------------------------------
#  compute the 1st deriv of array f assuming equally spaced points   5-pt formula
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft, rfft, dct
    
    N = f.size
    df = np.zeros_like(f)
    xx = 1./(12.*h)
    
    df[0] = xx*(-25.*f[0] +48.*f[1] -36.*f[2] +16.*f[3] -3.*f[4])
    df[1] = xx*(-25.*f[1] +48.*f[2] -36.*f[3] +16.*f[4] -3.*f[5])
    for i in range(2,N-2):
    	df[i]=xx*(f[i-2] -8.*f[i-1] +8.*f[i+1] -f[i+2])
    df[N-1] = xx*( 25.*f[N-1] -48.*f[N-2] +36.*f[N-3] -16.*f[N-4] +3.*f[N-5])
    df[N-2] = xx*( 25.*f[N-2] -48.*f[N-3] +36.*f[N-4] -16.*f[N-5] +3.*f[N-6])
    return df
