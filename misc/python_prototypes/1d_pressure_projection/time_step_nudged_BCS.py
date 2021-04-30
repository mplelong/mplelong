#!/opt/local/bin/python
"""
  Illustrate solving a problematic pressure poisson eqn for an open boundary
  problem when approach is limited for practical reasons to cos expansions of
  all dependent variables

    u_t = A*cos(2pi/L x) + nu u_xx  + BC_nudging(t) - p_x    u_initial = 1
       
    cos term models a body force that puts div into the flow but doesn't have a mean acceleration
    imagine that some large-scale external pressure gradient is slowly modifying the boundary vals
    how does our scheme manage this problem?
    
    u* = u_n + integral (body_force + diss + BC_nudging) dt    nudging term contains imposed boundary values
    div u* = d/dx u* 
    
    time step using simple Euler scheme with 
    u_BV = cos(2pi/T t)
    
    THIS IS A SIMPLIFICATION OF MY EXPLICIT BC INSERTION W/ BENDING NEAR EDGES
    TO JUST NUDGE THE BCS VIA FORCING        
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from KW_FUNCTIONS import dst,dst_filtered,compute_spectra,poisson_invert
#-----------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#  parameters for the example problem
#---------------------------------------------------------------------------------------
pi = np.pi
a = 0.0                 # [m]           left boundary position
b = 1.0                 # [m]           right boundary position
L = b-a          		# [m]			size of outer domain, soln periodic over L	
nx = 257                # [1]           number of discrete grid points, basis functions
M = 2*(nx-1)            # [1]           number of points in extended, open interval for periodic representation
dx = L/(nx-1.)          # [m]           grid spacing
dk = pi/L               # [1/m]         wavenumber interval
k_nyq = nx*pi/L			# [1/m]		    Nyquist wavenumber

                        

T = 1.0                 # [s] period of boundary forcing term
nt = 100
dt = T/nt               # [s]           time step   CFL u*dt < 1/2 dx   u_max < .5 dx/dt 
                        #               so, jumps in u at bdries for new information are much smaller that the max u
                        
nu = .000*(dx**2)/dt     # [m2/s]        viscosity
frac=0.05                # fraction of wavenumber space to filter when computing d/dx of u* and p   (must be > 0)

x = np.linspace(a,b,nx)             # [m] discrete grid points in [0.L]
x_ext = np.linspace(a,2*L-dx,M)     # [m] discrete grid points in [0,2L)
ustar = np.zeros_like(x)            # create some arrays for later use
ustar_ext = np.zeros_like(x_ext)

A = 1.0                               # [m/s2]  amplitude of body force
body_force = A*np.cos((2.*pi/L)*x)    #


#---------------------------------------------------------------------------------------
#  parameters for BC nudging window
#---------------------------------------------------------------------------------------
sigma = 5.*dx
W_left  = np.exp(-((x-a)/sigma)**4)    # near boundary nudging region  (4 is better than 2)
W_right = np.exp(-((x-b)/sigma)**4) 
tau = 3.*dt                                                 # nudging time scale


#---------------------------------------------------------------------------------------------
#  initial conditions  u(x) = 1.0
#---------------------------------------------------------------------------------------------
u_n = np.ones_like(x)


plot_step = 99
for j in range(nt):

	tn = j*dt
	tnp1 = (j+1)*dt
	u_BV = np.cos((2.*pi/T)*tnp1)      # BCs on u at t_n+1
	print('new BVAL: ',u_BV)
	
	#---------------------------------------------------------------------------------------------
	# use cos expansion of u* to compute its deriv d/dx u*
	#---------------------------------------------------------------------------------------------
	nderiv = 2                                   # take 2nd deriv
	flag = -1                                    # u_n even/cos expandable  (-1/1 for cos/sin expandable functions)
	u_xx,u_xx_ext = dst(u_n,nderiv,L,flag)       
	diss = nu*u_xx
	
	nudge = (-1./tau)*W_left*( u_n[:] - u_BV ) + (-1./tau)*W_right*( u_n[:] - u_BV )
	
	
	ustar = u_n + dt*( body_force + diss + nudge)    #  body force that generates divergence that needs elimination
	
	u_exact = u_BV*np.ones_like(x)


	#---------------------------------------------------------------------------------------------
	#  create the even extension of ustar [0,L] --> [0,2L)
	#---------------------------------------------------------------------------------------------
	ustar_ext = np.concatenate([ustar,  ustar[1:-1][::-1]])    # EVEN extension of u* to [0,2L)

	#---------------------------------------------------------------------------------------------
	# use cos expansion of u* to compute its deriv d/dx u*
	#---------------------------------------------------------------------------------------------
	nderiv = 1                                   # take 1st deriv
	flag = -1                                    # ustar even/cos expandable  (-1/1 for cos/sin expandable functions)
	#div,div_ext = dst(ustar,nderiv,L,flag)                # result is a NOT cos expandable function, div(u*)
	div,div_ext = dst_filtered(ustar,nderiv,L,flag,frac)
	rhs = div                                             # don't modify the computed div u* vector


	#---------------------------------------------------------------------------------------
	#  set up the figure for plotting u* and its even extension
	#---------------------------------------------------------------------------------------
	if( j==plot_step):
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
    
		title_string = r"$u^*$ obtained by integration with nudging term added"    
		axes.set_title(title_string,fontsize=12)
		axes.plot(x/L,ustar,'k',linewidth=1)
		axes.plot(x_ext/L,ustar_ext,'k',linewidth=1)
		axes.plot(x/L,ustar,'r.',markersize=6)
		axes.plot(x_ext/L,ustar_ext,'k.',markersize=4)
		axes.grid()
		
		#axes.plot(x/L,np.zeros_like(x),'k--',linewidth=0.5)
		
		plt.savefig('ustar_vecs.eps',dpi=300,bb_inches='tight')      # save plot file

	#---------------------------------------------------------------------------------------
	#  set up the figure for plotting the corresponding rhs vectors
	#---------------------------------------------------------------------------------------
	if( j==plot_step):
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
	if( j==plot_step):
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
	if( j==plot_step):
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
		axes.plot(k[1:-1:2]/dk,np.log10(np.abs(FHAT[1:-1:2]/k[1:-1:2]**2)),'ko',markersize=3)
		axes.plot(k[1:-1:2]/dk,np.log10(np.abs(FHAT[1:-1:2]/k[1:-1:2])),'ro',markersize=3)
		axes.grid()

		plt.savefig('p_spectrum.eps',dpi=300,bb_inches='tight')      # save plot file



	#---------------------------------------------------------------------------------------
	#  set up the figure for plotting the exact and computed solns
	#---------------------------------------------------------------------------------------
	if( j==plot_step):
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
    
		title_string = r"computed soln for pressure (k,r)"    
		axes.set_title(title_string,fontsize=12)
		axes.plot(x,p,'k',linewidth=1)
		axes.plot(x_ext,p_ext,'k',linewidth=1)
		axes.plot(x,p,'r.',markersize=6)
		axes.plot(x_ext,p_ext,'k.',markersize=4)

		#axes.plot(x,p_exact,'b',linewidth=2)
		#axes.plot(x_ext,rhs_ext/(2*pi),'g',linewidth=2)
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
	if( j==plot_step):
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
	#grad_p,grad_p_ext = dst(p,nderiv,L,flag)       # 2nd deriv of cos expandable function	
	grad_p,grad_p_ext = dst_filtered(p,nderiv,L,flag,frac)

	#---------------------------------------------------------------------------------------
	#  set up the figure for p_x
	#---------------------------------------------------------------------------------------
	if( j==plot_step):
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
    
		title_string = r"$\phi_x$  "    
		axes.set_title(title_string,fontsize=12)
		axes.plot(x/L,grad_p,'k',linewidth=1)
		axes.plot(x/L,grad_p,'r.',markersize=6)
		axes.plot(x_ext/L,grad_p_ext,'k',linewidth=1)
		axes.plot(x_ext/L,grad_p_ext,'k.',markersize=3)

		plt.savefig('pressure_gradient.eps',dpi=300,bb_inches='tight')      # save plot file


	#---------------------------------------------------------------------------------------
	#  set up the figure for pressure correction
	#---------------------------------------------------------------------------------------
	if( j==plot_step):
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
    
		title_string = r"$u(x)=u^* - \nabla (r)$  and $u_{\rm exact}$  (b)  "    
		axes.set_title(title_string,fontsize=12)
		#axes.plot(x/L,grad_p,'k',linewidth=1)
		#axes.plot(x/L,ustar,'k--',linewidth=1)
		axes.plot(x/L,ustar-grad_p,'r',linewidth=2)

		axes.plot(x/L,u_exact,'b*',markersize=1)
		#axes.plot(x/L,np.zeros_like(x),'k--',linewidth=0.5)
		axes.plot(x/L,1*np.ones_like(x),'k--',linewidth=0.5)
		axes.plot(x/L,np.cos(dt*2*pi/T)*np.ones_like(x),'k--',linewidth=0.5)
		axes.plot(x/L,.996*np.ones_like(x),'k--',linewidth=0.5)
		#axes.plot(x/L,-1*np.ones_like(x),'k--',linewidth=0.5)

		plt.savefig('pressure_correction.eps',dpi=300,bb_inches='tight')      # save plot file
		
		error = np.mean(np.abs((ustar-grad_p)-u_exact))
		print(error)
		exit()

	#---------------------------------------------------------------------------------------
	#  do the pressure projection
	#---------------------------------------------------------------------------------------
	u_n = ustar - grad_p
	
exit()



