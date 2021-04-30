#!/opt/local/bin/python
"""
    1d wave equation example   u_t + cu_x = 0   periodic in [0,2pi], c specified, inviscid
    large scale solution: u = g(x-ct)    ; let g(xsi) = sin(xsi) for example
    soln satisfies eqn and periodic BCs at x=0, x=2pi
    
    try to recover solution within subdomain, 
    using the pde and "BCs" at 2 arbitrary locations  x0 and x1
    try to use the minimal diffusion necessary to maintain stability
    
    soln is neither periodic nor even/odd symmetric over this subdomain
    but we will do the numerical integration using cosine expansions and
    boundary information supplied externally, i.e. extracted from the
    global problem
    
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from KW_FUNCTIONS import dst,dst_filtered,diffusion_2p,compute_spectra,spectral_filter
#-----------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
#  For the record, a parameter set that gives really nice solns:
#    N=1025, 
#    dt = 0.06*dx/c 
#    nu_star = 0.004 * (dx^2p)/dt
#    p=1,2,3   (p=3 also worked w/ same parameters and N=33,65,129,257,513) stable but inaccurate @33,65
#    beta=sigma=0.001dx  Neumann correction False
#    2p spectral filtering highest 10% of wavenumber space 
#    alpha=pi/3. for exponential function g(xsi)
#
#    N=2049 worked with above except dt=.04*dx/c
#    p=4 N=513 .06dx/c 30% filter   nu_star = 0.004 * (dx^2p)/dt, 0.001 * (dx^2p)/dt (bigger/smaller unstable)
#---------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#  parameters for the example problem
#---------------------------------------------------------------------------------------
pi = np.pi
L_outer = 2.0*pi		# [m]			size of outer domain, soln periodic over L	
x0 = 0.9*pi				# [m]			left edge of inner subdomain
x1 = 1.1*pi				# [m]			right edge of inner subdomain
L = x1 - x0             # [m]			inner domain length   L << L_outer

g_type = 'pulse'   		# 				'pulse' or 'plane_wave'
alpha = pi/3.          	# [m]       	spatial scale for pulse wave function, outer solution
kp = 2.0				# [1/m]			spatial wavenumber of plane wave outer solution
c = 1.0					# [m/s]			prescribed wave speed
T_traverse = L_outer/c	# [s]			time for traversal of outer domain
T_f = 0.50*T_traverse	# [s]			integration length   

N_outer = 257								#			outer soln and thus BCs known analytically
dx_outer = L_outer/(N_outer-1.)				# [m]		interior grid spacing, closed domain
x_outer = np.arange(N_outer)*dx_outer		# [m]		discrete interior x gridpoints, including endpoints

N = 513  						# []		number of gridpoints in interior subdomain 
dx = L/(N-1.)					# [m]		interior grid spacing, closed domain
x = x0 + np.arange(N)*dx		# [m]		discrete interior x gridpoints, including endpoints

k_nyq = N*pi/L					# [1/m]		Nyquist wavenumber for interior subdomain
dt = 0.20*dx/c					# [s]		discrete time step
nt = np.int(T_f/dt)				# []		integer number of time steps to take
t = np.arange(nt)*dt			# [s]		array of discrete time values

p=4                     		# [1]       2p is order of diffusion operator
nu_star = (.00014)*(dx**(2*p))/dt   # [m**2p/s]	diffusion coefficient   (p=1;3.25e-5 m2/s)
frac = 0.125					# 			fraction of high wavenumbers to filter from derivatives/computed soln  3,.15

beta = 0.00001*dx				# [m]		spatial scale for local boundary value corrections beta = O(dx)  1.5
sigma = beta					# [m]		spatial decay scale for perturbations needed for Neumann conditions
neumann_correction = False      	# whether or not to add perturbation to ensure no normal derivs at bdries

print "p,nu_*, dt, dx:  ",p,nu_star,dt,dx
print "stability ratios: adv,diff  ", c*dt/dx , nu_star*dt/dx**(2*p)
if(nu_star*dt/dx**(2*p) > .25):
	print "too much diffusion: nu_star*dt/dx**2p > .25  ", nu_star*dt/dx**(2*p)
	exit()
if(nu_star*dt/dx**(2*p) < 1.e-19):
	print "too little diffusion: nu_star*dt/dx**(2*p) < .0001  ", nu_star*dt/dx**(2*p)
	exit()
if(c*dt/dx > 0.22):    # NICE TO HAVE THIS VALUE EXPERIMENTALLY!
	print "advective stability problem"
	exit()

	
#---------------------------------------------------------------------------------------
#  initial conditions at discrete gridpoints in inner domain, t=0
#  set u = g(x-ct)    2 sample functions for g
#---------------------------------------------------------------------------------------
A = 1.0
if( g_type=='plane_wave' ):
	xsi = kp*(x - c*t[0])
	u = A*np.sin( xsi )     			# g = sin(xsi)
if( g_type=='pulse' ):
	xsi = x - c*t[0]
	u = A*np.exp(-(xsi/alpha)**2)		# g = exp(-(xsi/alpha)**2)

#---------------------------------------------------------------------------------------
# generate BC values at x0,x1 for each discrete time step
#   (this mimics saving boundary values from an outer lower resolution simulation)
#---------------------------------------------------------------------------------------
if( g_type=='plane_wave' ):
	xsi = kp*(x0 - c*t)
	U0 = A*np.sin( xsi )
if( g_type=='pulse' ):
	xsi = x0 -c*t  					# g = sin(xsi)
	U0 = A*np.exp(-(xsi/alpha)**2)	# g = exp(-(xsi/alpha)**2)



if( g_type=='plane_wave' ):
	xsi = kp*(x1 - c*t)
	U1 = A*np.sin( xsi )
if( g_type=='pulse' ):
	xsi = x1 -c*t   				# g = sin(xsi)
	U1 = A*np.exp(-(xsi/alpha)**2)	# g = exp(-(xsi/alpha)**2)



#---------------------------------------------------------------------------------------
#  set up the figure for plotting the solution
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(16,8),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x$ [m]',fontsize=14)
axes.set_ylabel(r'$u$ [m/s]',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     	# Get rid of top axis line
axes.spines['right'].set_visible(False) 
    
title_string = "1-dimensional wave equation"    
axes.set_title(title_string,fontsize=12)


#---------------------------------------------------------------------------------------
#  loop over discrete time steps
#---------------------------------------------------------------------------------------
for i in range(nt):						# range(nt)
		
	#--------------------------------------------------------
	#  evaluate outer solution  (exact wave eqn soln)
	#--------------------------------------------------------	
	if( g_type=='plane_wave' ):
		xsi = kp*(x_outer - c*t[i])
		u_outer = A*np.sin( xsi )    			# g = sin(xsi)
	if( g_type=='pulse' ):
		xsi = x_outer - c*t[i]
		u_outer = A*np.exp(-(xsi/alpha)**2)		# g = exp(-(xsi/alpha)**2)
	
	#--------------------------------------------------------
	#  impose inhomogeneous Dirichlet conditions at edges
	#   use local corrections over spatial scale beta
	#   BY EXPERIMENTATION: beta can --> 0 
	#--------------------------------------------------------
	#u_d = u + (U0[i]-u[0])*np.exp(-(x-x0)/beta) + (U1[i]-u[-1])*np.exp(-(x1-x)/beta)
	u_d = u
	u_d[0] = U0[i]
	u_d[-1] = U1[i]
	
	#--------------------------------------------------------
	#  impose homogeneous Neumann conditions at edges
	#  by adding small global perturbations
	#    BY EXPERIMENTATION: THIS STEP IS NOT NECESSARY
	#--------------------------------------------------------
	
	#------------------------------------------------------
	# (a)  fix (or select nearby value stochastically) k,
	#      the wavenumber of the perturbation field
	#      k should be O(k_nyq)
	#------------------------------------------------------
	k = .975*k_nyq    # .975     NOT ACTUALLY USED
			
	#------------------------------------------------------
	# (b)  estimate derivs of the function u_d at x0,x1
	#      use these and k to compute amplitudes of corrections
	#      needed to enforce homogeneous Neumann at endpoints
	#------------------------------------------------------
	deriv_0 = (u_d[1]-u_d[0])/dx		# simple forward difference
	deriv_1 = (u_d[-1]-u_d[-2])/dx		# simple backward difference
	A0 = -deriv_0/k
	A1 = -deriv_1/k
		
	#------------------------------------------------------
	# (c)  add small amplitude corrections to u_d, 
	#	   localize near boundary using scale sigma
	#      to get u_dn; a perturbed version of the desired
	#      solution satisfying the D BCs but expandable in cos
	#------------------------------------------------------
	if(neumann_correction):
		u_dn = u_d + A0*np.sin(k*(x-x0))*np.exp(-((x-x0)/sigma)**2) + A1*np.sin(k*(x-x1))*np.exp(-((x-x1)/sigma)**2)
	else:
		u_dn = u_d   # JUST USE U W/ BDRY VALUES REPLACED
	
	
	#--------------------------------------------------------
	#  now the time stepping can be undertaken with spatial
	#  derivatives computed spectrally
	#--------------------------------------------------------
	n=1   											# 1st deriv
	flag = -1										# function is even, i.e. cos expandable
	dudx = dst_filtered(u_dn,n,L,flag,frac)			# compute spatial deriv w/ cos transform			
	diff_term = diffusion_2p(u_dn,p,L,flag,nu_star,frac)    # (-1)**(p-1) * nu_star * 2p deriv of u_dn, deriv is filtered
	
	rhs_i = -c*dudx + diff_term		# compute current rhs
	
	if(i==0):
		u_new = u_dn + dt*rhs_i															# 1st step Euler
	elif(i==1):
		u_new = u_dn + (dt/2.)*(3.*rhs_i - rhs_im1)										# 2nd step AB2
	elif(i==2):
		u_new = u_dn + (dt/12.)*(23.*rhs_i - 16.*rhs_im1 + 5.*rhs_im2)					# AB3
	elif(i>=3):
		u_new = u_dn + (dt/24.)*(55.*rhs_i - 59.*rhs_im1 + 37.*rhs_im2 - 9.*rhs_im3)	# AB4
	
	#---------	
	# toggle
	#---------
	rhs_im1 = rhs_i	
  	rhs_im2 = rhs_im1
  	rhs_im3 = rhs_im2
  	
  
  
	#---------------------------------------------------------
	# update the u_i to u_i+1
	# by saving spectrally filtered time stepped solution
	#---------------------------------------------------------
	flag = -1               # u is even/cos
	u = spectral_filter(u_new,frac,L,flag,p)   #  spectrally filter and go on to next time step (i loop)
	#u = u_computed	
	
	
	if( np.mod(i,np.int(nt/16))==0 ):
		print "plotting step ",i,t[i]
		print u[0],u[-1]
		#----------------------------------------------------------------
		#  plot 
		#   u_outer    		exact soln over outer/global domain
		#	u     	 		new soln after time step
		#----------------------------------------------------------------
		plt.plot(x,u,'b',linewidth=1.25,label=r'$u_{computed}$')
		plt.plot(x_outer,u_outer,'r',linewidth=0.5,label=r'$u_{outer}$')

		

#------------------------------------------------------
#	add subdomain boundaries to figure and save as eps
#------------------------------------------------------
if( g_type=='plane_wave' ): y0=-A
if( g_type=='pulse' ): 		y0=0
plt.plot([x0,x0],[y0,A],'b--',linewidth=1)
plt.plot([x1,x1],[y0,A],'b--',linewidth=1)
plt.savefig('test.eps',dpi=300,bb_inches='tight')      # save plot file



#---------------------------------------------------------------------------------------
#  set up the figure for plotting spectra
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(16,8),dpi=300)              # fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.10, 0.85, 0.85])    # lower, bottom, width, height in ( 0 to 1)


[k,E,D] = compute_spectra(u,flag,L,nu_star,p)
plt.subplot(121)
plt.loglog(k,E)
plt.loglog([(1-frac)*k[-1],(1-frac)*k[-1] ],[1.e-11,1.e-7],'r')
plt.loglog([.95*k[-1],.95*k[-1] ],[1.e-11,1.e-7],'r')
plt.subplot(122)
plt.loglog(k,1/D)
plt.savefig('spectra.eps',dpi=300,bb_inches='tight')      # save plot file

exit()



