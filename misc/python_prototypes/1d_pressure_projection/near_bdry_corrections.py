#!/opt/local/bin/python
"""
  Illustrate solving a problematic pressure poisson eqn for an open boundary
  problem when approach is limited for practical reasons to cos expansions of
  all dependent variables

     p_xx = f(x)                  f =   x^3 - L^3/4   
                        
    BCS:  p_x(0)=p_x(L) = 0     i.e. BCs seemingly consistent w/ p being expandable in cos series

     exact solution:  p(x) = (1/20) x^5 - (L^3/8) x^2
     
            	==>   p_x  = (1/4) x^4 - (L^3/4) x   p_x(0) = 0, p_x(L) = 0

     this script will document the approach and quantify how good an approximation
     to the solution we are able to abtain
    
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
def poisson_invert(f,L):
#  solve p_xx = f(x)   p_x(0)=p_x(L)=0, equally spaced closed interval grid p,f, even
#  flag =-1  ==> assume that f is expandable in a cos series, 
#                f(-x)=f(x) near x=0,L  f(0)=f(L)  explicitly stored
#                work with explicit even extension in Fourier space for [0,L)
#-----------------------------------------------------------------------------------------
    import numpy as np
    from scipy.fftpack import fft, ifft
    
    N = f.size
    dk = np.pi/L
    M = (N-1)*2                                              # extended array length
    k = dk * np.array( range(M/2+1) + range(-M/2+1,0,1) )
    F = np.concatenate([f,  f[1:-1][::-1]])                  # F is even extension of data
    FHAT = fft(F)                                            # FT of even extended periodic data vector
    phat = np.zeros_like(FHAT)
    
    for i in range(1,M):                                     # skip 0 wavenumber location
    	phat[i] = -(1./k[i]**2) * FHAT[i]                    # integrate twice wrt x in wavenumber space
        
	p = ifft(phat)[0:M/2+1].real                             # inverse transform, keep cos series
    return p
#-----------------------------------------------------------------------------------------






#---------------------------------------------------------------------------------------
#  parameters for the example problem
#---------------------------------------------------------------------------------------
pi = np.pi
a = 0.0                 #  left boundary position
b = 1.0                 #  right boundary position
L = b-a          		#  size of outer domain, soln periodic over L	
nx = 257                #  number of discrete grid points, basis functions

x = np.linspace(a,b,nx) # [m]           discrete grid points in [0.L]

#------------------------------------------------------------------------
#  this problem is treated as even symmetric but the even extension has 
#  a discontinuous derivative at the end points
#f = x**3 - (L**3)/4.                          # rhs  compatable w/ homog. Neumann
#phi_exact = (1./20.)*x**5 - ((L**3)/8.)*x**2


#------------------------------------------------------------------------
#  this problem is truly even symmetric
#f = np.cos(x*pi/L)                           #  compatable
#phi_exact = -(L/pi)**2 * np.cos(x*pi/L)


#------------------------------------------------------------------------
#  this problem is treated as even symmetric but the even extension has 
#  a discontinuous derivative at the end points
k = 2.*pi/L ; phase = pi/4. ; arg = k*x + phase ; C = (1./k)*np.cos(phase)
f = np.sin(arg)                        #  compatable
phi_exact = -(1/k)**2 * np.sin(arg) + C*x

x0 = 0.0 ; arg = k*x0 + phase
f0 = np.sin(arg)             # f(0)
fp = k*np.cos(arg)           # f'(0)
fpp = -k*k*np.sin(arg)       # f"(0)
fppp = -k*k*k*np.cos(arg)    # f"'(0)   etc



#---------------------------------------------------------------------------------------
#  invert the Poisson eqn using even cos expansions
#---------------------------------------------------------------------------------------
phi = poisson_invert(f,L)
phi = (phi-phi[0]) + phi_exact[0]


#---------------------------------------------------------------------------------------
#  compute the error and even and odd terms in Taylor series expansion
#---------------------------------------------------------------------------------------

error = phi_exact - phi

x0 = 0. ; C = 0.
series = C*x + (f0/2.)*(x-x0)**2 + (fp/6.)*(x-x0)**3 + (fpp/12.)*(x-x0)**4 + (fppp/120.)*(x-x0)**5
even_terms = (f0/2.)*(x-x0)**2 + (fpp/12.)*(x-x0)**4
odd_terms = C*x + (fp/6.)*(x-x0)**3 + (fppp/120.)*(x-x0)**5 


#---------------------------------------------------------------------------------------
#  set up the figure for plotting phi_even and phi_exact
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r'',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"$\phi$ and $\phi_{\rm even}$"    
axes.set_title(title_string,fontsize=12)

npts=5


axes.plot(x,error,'k',linewidth=1)
#axes.plot(x[0:npts],even_terms[0:npts],'b*',linewidth=1)
axes.plot(x,phi_exact,'r',x,phi,'k',linewidth=1)
axes.grid()

plt.savefig('test.eps',dpi=300,bb_inches='tight')      # save plot file



exit()



