#!/opt/local/bin/python
"""
  Illustrate solving a problematic pressure poisson eqn for an open boundary
  problem when approach is limited for practical reasons to cos expansions of
  all dependent variables

    u_t = cos(2pi/L x) - p_x 
       
    cos term models a body force that puts div into the flow but doesn't have a mean acceleration
    imagine that some large-scale expernal pressure gradient is slowly modifying the boundary vals
    how does our scheme manage this problem?
    
    u* = u_n + dt*cos(2pi/L x) in the interior  BV(t) at x=0,L   bdry vals BV obtained externally
    div u* = d/dx u* = -(dt*2pi/L)*sin(2pi/L x)
    
    
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from KW_FUNCTIONS import dst,compute_spectra,poisson_invert
#-----------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#  parameters for the example problem
#---------------------------------------------------------------------------------------
pi = np.pi
a = 0.0                 # [m]           left boundary position
b = 1.0                 # [m]           right boundary position
L = b-a          		# [m]			size of outer domain, soln periodic over L	
nx = 65                 # [1]           number of discrete grid points, basis functions
M = 2*(nx-1)            # [1]           number of points in extended, open interval for periodic representation
dx = L/(nx-1.)          # [m]           grid spacing
gamma = 2.*dx    		# [m]           decay length for "bending" near boundaries  (require gamma<<L)
dk = pi/L               # [1/m]         wavenumber interval
k_nyq = nx*pi/L			# [1/m]		    Nyquist wavenumber
dt = 0.1                # [s]           time step   CFL u*dt < 1/2 dx   u_max < .5 dx/dt = .5*L/((nx-1)*.1) ~ 0.15 m/s
                        #               so, jumps in u at bdries for new information are much smaller that the max u

x = np.linspace(a,b,nx)             # [m] discrete grid points in [0.L]
x_ext = np.linspace(a,2*L-dx,M)     # [m] discrete grid points in [0,2L)
ustar = np.zeros_like(x)
ustar_ext = np.zeros_like(x_ext)

u_initial = np.ones_like(x)         # u_initial, i.e. u at t_n
A = 0.01                            # amplitude of "body force"

#---------------------------------------------------------------------------------------------
#  define ustar(x) as if it came from a body force that creates some divergence but no mean accel
#  and then impose the BCs for u at time (n+1) onto u* itself, i.e. small jumps at boundaries
#---------------------------------------------------------------------------------------------
ustar = u_initial + A*np.cos((2.*pi/L)*x)     # w/o updated boundary information


C = 0.95*u_initial[0]                         # jumps in u at boundaries imposed from external information
ustar[0]=C ; ustar[-1]=C                      # BVALS
u_exact = C*np.ones_like(x)                   # divergence free solution that satisfies BCs on u


#---------------------------------------------------------------------------------------
#  bend the u* function near the edges to give it zero derivatives there...
#  call the result rhs (or, test no bending)
#---------------------------------------------------------------------------------------
for i in range(nx):
	if( x[i] <= L/2. ):
		g = 1.0 - np.exp(-((x[i]-a)/gamma)**2)
		ustar[i] = (ustar[i]-ustar[0])*g + ustar[0]
	else:
		g = 1.0 - np.exp(-((x[i]-b)/gamma)**2)
		ustar[i] = (ustar[i]-ustar[-1])*g + ustar[-1]


#---------------------------------------------------------------------------------------------
#  create the even extension of ustar [0,L] --> [0,2L)
#---------------------------------------------------------------------------------------------
ustar_ext = np.concatenate([ustar,  ustar[1:-1][::-1]])    # EVEN extension of u* to [0,2L)

#---------------------------------------------------------------------------------------------
# use cos expansion of u* to compute its deriv d/dx u*
#---------------------------------------------------------------------------------------------
nderiv = 1                                   # take 1st deriv
flag = -1                                    # ustar even/cos expandable  (-1/1 for cos/sin expandable functions)
div,div_ext = dst(ustar,nderiv,L,flag)       # result is a NOT cos expandable function, div(u*)
rhs = div                                    # don't modify the computed div u* vector


#---------------------------------------------------------------------------------------
#  set up the figure for plotting u* and its even extension
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
    
title_string = r"$u^*$ with bdry values imposed and smoothed, $\gamma=2\Delta x$"    
axes.set_title(title_string,fontsize=12)
axes.plot(x/L,ustar,'k',linewidth=1)
axes.plot(x_ext/L,ustar_ext,'k',linewidth=1)
axes.plot(x/L,ustar,'r.',markersize=6)
axes.plot(x_ext/L,ustar_ext,'k.',markersize=4)
axes.grid()

plt.savefig('ustar_vecs.eps',dpi=300,bb_inches='tight')      # save plot file

#---------------------------------------------------------------------------------------
#  set up the figure for plotting the corresponding rhs vectors
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
    
title_string = r"$u^*_x$"    
axes.set_title(title_string,fontsize=12)
axes.plot(x/L,div,'k',linewidth=1)
axes.plot(x_ext/L,div_ext,'k',linewidth=1)
axes.plot(x/L,div,'r.',markersize=6)
axes.plot(x_ext/L,div_ext,'k.',markersize=4)
axes.grid()

plt.savefig('rhs_vecs.eps',dpi=300,bb_inches='tight')      # save plot file

#---------------------------------------------------------------------------------------
#  invert the modified rhs to generate an approximate solution
#---------------------------------------------------------------------------------------
p,p_ext,rhs_ext,FHAT,k = poisson_invert(rhs,L)


#---------------------------------------------------------------------------------------
#  set up the figure for plotting Fourier spectrum of rhs in [0,2L) w/ discontinuous derivative
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$k/\Delta k$ [m]',fontsize=14)
axes.set_ylabel(r'$\log_{10} |\hat f|$',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"FFT of even-extended rhs of p eqn w/ discontinuous derivative (every other coeff = 0)"    
axes.set_title(title_string,fontsize=10)
axes.plot(k[1:-1:2]/dk,np.log10(np.abs(FHAT[1:-1:2])),'ko',markersize=3)
axes.grid()

plt.savefig('rhs_spectrum.eps',dpi=300,bb_inches='tight')      # save plot file


#---------------------------------------------------------------------------------------
#  set up the figure for plotting Fourier spectrum of p_ext in [0,2L) w/ discontinuous derivative
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$k/\Delta k$ [m]',fontsize=14)
axes.set_ylabel(r'$\log_{10} |\hat f/k^2|$',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"${\hat p}$ (black)            $k{\hat p} \sim {\hat u}$  (red),   (every other coeff = 0)"    
axes.set_title(title_string,fontsize=12)
axes.plot(k[1:-1:2]/dk,np.log10(np.abs(FHAT[1:-1:2])/k[1:-1:2]**2),'ko',markersize=3)
axes.plot(k[1:-1:2]/dk,np.log10(np.abs(FHAT[1:-1:2])/np.abs(k[1:-1:2])),'ro',markersize=3)
axes.grid()

plt.savefig('p_spectrum.eps',dpi=300,bb_inches='tight')      # save plot file



#---------------------------------------------------------------------------------------
#  set up the figure for plotting the exact and computed solns
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x$ [m]',fontsize=14)
axes.set_ylabel(r'',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     	# Get rid of top axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"computed soln (k,r),  scaled, extended rhs (g)"    
axes.set_title(title_string,fontsize=12)
axes.plot(x,p,'k',linewidth=1)
axes.plot(x_ext,p_ext,'k',linewidth=1)
axes.plot(x,p,'r.',markersize=6)
axes.plot(x_ext,p_ext,'k.',markersize=4)

#axes.plot(x,p_exact,'b',linewidth=2)
axes.plot(x_ext,rhs_ext/(2*pi),'g',linewidth=2)
axes.grid()

plt.savefig('p_soln.eps',dpi=300,bb_inches='tight')      # save plot file


#---------------------------------------------------------------------------------------
#  compute d_xx of p in [0,L] using cos expansion
#---------------------------------------------------------------------------------------
nderiv = 2; flag = -1
pxx,pxx_ext = dst(p,nderiv,L,flag)       # 2nd deriv of cos expandable function
max = np.max( np.abs(div) )

#---------------------------------------------------------------------------------------
#  set up the figure for Residual = pxx - div    (should be 0 for exact solution p)
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r'',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     	# Get rid of top axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"pointwise residual $(p_{xx} - u^*_x)/{\rm max}(u^*_x)$ "    
axes.set_title(title_string,fontsize=12)
axes.plot(x/L,(pxx-div)/max,'k',linewidth=1)
axes.plot(x/L,(pxx-rhs)/max,'r.',Markersize=0.5)

plt.savefig('residual.eps',dpi=300,bb_inches='tight')      # save plot file



#---------------------------------------------------------------------------------------
#  compute grad p = d/dx of p in [0,L] using cos expansion
#---------------------------------------------------------------------------------------
nderiv = 1; flag = -1
grad_p,grad_p_ext = dst(p,nderiv,L,flag)       # 2nd deriv of cos expandable function


#---------------------------------------------------------------------------------------
#  set up the figure for Residual = pxx - div    (should be 0 for exact solution p)
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r'',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     	# Get rid of top axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"$\nabla p = p_x~(k),~u^*~(k--)$ and correction $u(x)=u^* - \nabla (r)$  "    
axes.set_title(title_string,fontsize=12)
axes.plot(x/L,grad_p,'k',linewidth=1)
axes.plot(x/L,ustar,'k--',linewidth=1)
axes.plot(x/L,ustar-grad_p,'r',linewidth=1)

axes.plot(x/L,u_exact,'b*',markersize=3)
axes.plot(x/L,u_initial,'g--',linewidth=1)

plt.savefig('pressure_correction.eps',dpi=300,bb_inches='tight')      # save plot file




exit()



