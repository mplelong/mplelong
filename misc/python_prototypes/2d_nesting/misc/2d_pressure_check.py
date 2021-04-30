#!/usr/local/bin/python
"""
  Everything SEEMS clean, step by step but the boundaries are not handled well
  
  
  Illustrate solving for 2d flow in a nested domain with open boundaries.
  BCs are obtained from analytic expressions for an internal wave mode propagating 
  through an arbitrary interior x-z subdomain of interest. The full domain has a 
  solid bottom and rigid lid w/ free-slip, no normal flow bcs on the outer soln.
  
  Numerics:
  (1) Bernoulli/cosine method for derivatives, no symmetries assumed.
  (2) cosine-transform based fast solver for Poisson w/ homogeneous Neumann BCs
  (3) BCs on normal flow components embedded into source equation for pressure
  (4) BCs on normal derivatives of tangential velocity components also embedded
  
  ===> still need scheme for high, even derivatives (orders 4,6,8 say)

  2D internal wave mode:
  	U = A cos(m*Z)*cos(k*X-omega*t) 
  	V = A*(f/omega)*cos(m*Z)*sin(k*X-omega*t) 
  	W = A*(k/m)*sin(m*Z)*sin(k*X-omega*t) 
  	B = -A*(k/m)(N^2/omega)*sin(m*Z)*cos(k*X-omega*t) 
  	with omega^2 = (k^2*N^2 + m^2*f^2)/(k^2 + m^2) 
  	m = (pi/H)*{1,2,3...}    k = (2pi/L)*{1,2,3...} 
  	NB  B=W=0 at z=0,H    flow is periodic in x with periodicity L
  	
  Governing equations in 2d:  
       w/ 1/rho_0 absorbed into definition of p
       b = -g/rho_0 * rho'   rho = rho_bar(z) + rho'
       N^2 = -(g/rho_0)*d/dz(rho_bar) = constant > 0
       
  	u_t =  fv - (u u_x + w u_z) + nu Laplacian(u) - p_x          
  	v_t = -fu - (u v_x + w v_z) + nu Laplacian(v)
  	w_t =   b - (u w_x + w w_z) + nu Laplacian(w) - p_z            
  	b_t = - N^2 w - (u b_x + w b_z) + kappa Laplacian(b)           
  	u_x + w_z = 0
  	  
  BCs:
   [u,v,w,b](t) specified by evaluating [U,V,W,B] on boundary of nested domain
   
  Nested domain: 
  	x in [0,Lx]  z in [0,Lz]  offset from outer domain by x0 and z0
  	==> X = x0 + x  and  Z = z0 + z
  	
  2D array storage convention:
   u[x,z] , number of gridpoints on closed interior doman = nx and nz
   
  variable id = [0,1,2,3,4] for [u,v,w,b,p] respectively
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from solver_utilities import step, exp, create_inversion_matrix, poisson_invert_2d
from solver_utilities import prepare_BernoulliCosine, evaluate_basis_functions
from solver_utilities import create_2d_variables, create_1d_arrays
from solver_utilities import time_step, divergence, grad, grad2, endpoint_deriv

from netcdf_stuff import write_netcdf_2d


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#  parameters for the environment and the outer solution which is used to
#  supply initial and time dependent boundary conditions
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
pi = np.pi
H = 3000.                 # [m]           full water depth
L = 1.5e5                 # [m]           horizontal scale = 150 km
N = 2.0e-3                # [1/s]         buoyancy frequency 
BVPER = (2.*pi/N)/3600.   # [hours]       buoyancy period
f = 1.e-4                 # [1/s]         Coriolis parameter
IPER = (2.*pi/f)/3600.    # [hours]       inertial period
N2 = N*N; f2 = f*f        # [1/s2]        squared frequencies
nu = 0.e-5                # [m2/s]        viscosity
kappa = nu                # [m2/s]        diffusivity

k_iw = 2.*pi/L            # [1/m]         horizontal wavenumber of internal wave
m_iw = 1.*pi/H            # [1/m]         vertical wavenumber of internal wave

omega2 = (k_iw**2 * N2 + m_iw**2 * f2)/(k_iw**2 + m_iw**2)

omega = np.sqrt(omega2)   # [1/s]         internal wave frequency
IWPER = (2.*pi/omega)     # [s]
A = 0.01                  # [m/s]         amplitude of U
phase = pi/4.             # [1]           initial phase shift of wave mode


def parent_soln(X,Z,t,id):
	# function to evaluate parent soln at X,Z,t  ; 
	# here an IW mode with wave parameters available globally
	# id = [1,2,3,4,5] for [U,V,W,B,P] respectively
	argz = m_iw*Z
	argx = k_iw*X - omega*t + phase
	if(id==0): ans = A*np.cos(argz)*np.cos(argx) # U
	if(id==1): ans = A*(f/omega)*np.cos(argz)*np.sin(argx) # V 
	if(id==2): ans = A*(k_iw/m_iw)*np.sin(argz)*np.sin(argx) # W 
	if(id==3): ans = -A*(k_iw/m_iw)*(N2/omega)*np.sin(argz)*np.cos(argx) # B
	if(id==4): ans = -A*(1./(k_iw*omega))*(f2-omega2)*np.cos(argz)*np.cos(argx) # P	
	return ans

def parent_derivs(X,Z,t,dir):
	# function to evaluate certain derivs of parent soln at X,Z,t  ; 
	# here an IW mode with wave parameters available globally
	# dir = 'x','y' or 'z'
	# used for normal derivs of tangential vels at boundaries
	argz = m_iw*Z
	argx = k_iw*X - omega*t + phase
	
	if(dir=='x'):    # return V_x,W_x  (in this order)
		ans1 = k_iw*A*(f/omega)*np.cos(argz)*np.cos(argx)   # V_x
		ans2 = k_iw*A*(k_iw/m_iw)*np.sin(argz)*np.cos(argx) # W_x
	if(dir=='y'):    # return U_y,W_y 
		ans1 = 0.
		ans2 = 0.
	if(dir=='z'):    # return U_z,V_z
		ans1 = -m_iw*A*np.sin(argz)*np.cos(argx)             # U_z
		ans2 = -m_iw*A*(f/omega)*np.sin(argz)*np.sin(argx)   # V_z		
	return ans1,ans2

	

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#  parameters to define and control the nested simulation
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
Lx = L/5.                # [m]      size of nested domain in x  10 km
Lz = H/5.                # [m]      size of nested domain in z
x0 = 0.50*L				 # [m]      horizontal offset for origin of nested domain
z0 = 0.60*H              # [m]      vertical offset for origin of nested domain

# x nested grid arrays	
nx = 257                            # [1] number of discrete grid points, basis functions
dx = Lx/(nx-1.)                     # [m] grid spacing
x = np.linspace(0,Lx,nx)            # [m] discrete grid points in [0,Lx]
gamma_x = 4.*dx                     # [m] width of near boundary transition region   8dx


# z nested grid arrays	
nz = 257                            # [1] number of discrete grid points, basis functions
dz = Lz/(nz-1.)                     # [m] grid spacing
z = np.linspace(0,Lz,nz)            # [m] discrete grid points in [0,Lz]
gamma_z = 4.*dz                     # [m] width of near boundary transition region

#------------------------------------------------------------
# parameters controlling the discrete time integration
#------------------------------------------------------------
nt = 1500                           # [1] time steps per wave period
dt = IWPER/nt                       # [s] time step 
AB_order = 4                        # order of Adams Bashforth time stepping... (2 <= AB_order <= 4)
nsteps = 5
inc_write_netcdf = 1

#------------------------------------------------------------
# numerical differentiation parameters
#------------------------------------------------------------
Q = 7        #  ==> (Q+1)/2 = 4 terms in Bernoulli expansions
frac = 0.    # Fourier wavenumber filtering parameter


print( "dx  in m             ",dx)
print( "dz  in m             ",dz)
print( "dt  in seconds       ",dt)
print( "cfl x  U dt/dx       ",A*dt/dx )
print( "cfl z  W dt/dz       ",A*(k_iw/m_iw)*dt/dz )
print( "IWPER  in hours      ",IWPER/3600.) 
print( "IPER   in hours:     ",IPER)
print( "BVPER  in mins:      ",BVPER*60)
print( "omega/f              ",omega/f)


#-------------------------------------------------------------------------------------
#     PRELIMINARY TASKS
#-------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
# create the (extended) inversion matrix  -1/k2  for Fourier pressure solve 
# (actually applied in 2d wavenumber space to even extended div ustar)
#--------------------------------------------------------------------------
inversion_matrix = create_inversion_matrix(nx,nz,Lx,Lz)

#--------------------------------------------------------------------------
# evaluate the basis functions for the S_0 and S_L (U_n/B_n+1) series
# and for their derivatives 
#   x_basis_functions=[S_0,S_L,S_0_x,S_L_x] similar for z
#--------------------------------------------------------------------------
x_basis_functions = evaluate_basis_functions(x,Q)
z_basis_functions = evaluate_basis_functions(z,Q)
basis_functions = [x_basis_functions,z_basis_functions]

#--------------------------------------------------------------------------
# create and factor the matrices for solving for the expansion coeffs
# in the Bernoulli/Cosine differentiation scheme
#--------------------------------------------------------------------------
LU_x = prepare_BernoulliCosine(x,Q)
LU_z = prepare_BernoulliCosine(z,Q)


#---------------------------------------------------------------------------------------
#  construct some useful arrays
#---------------------------------------------------------------------------------------
#  step_x is an array that is 1 at the x boundaries decays to interior w/ 0 slope at ends
#  exp_x is an array that decays from E/W boundaries with slope at ends \pm 1/gamma_x
step_x_e,step_x_w  = step(x,gamma_x)
exp_x_e ,exp_x_w   =  exp(x,gamma_x)

# similar for z but with subscripts b & t for bottom and top
step_z_b,step_z_t = step(z,gamma_z)
exp_z_b ,exp_z_t  = exp(z,gamma_z)

# initialize some 2d arrays with zeros
[u,v,w,b,ustar,vstar,wstar,g2u,g2v,g2w,g2b] = create_2d_variables(nx,nz)

# initialize 4 time levels of the 2d rhs for each of the 4 variables u, v, w and b
RHS = np.zeros((nx,nz,4,4),float)     #  [nx,nz,var_id,4] last index for 4 time levels in AB4

#initialize the cyclic array for defining previous time levels          
idx = [0,1,2,3]   # i.e. [MM0,MM1,MM2,MM3], will be toggled in the time_step routine

# 1d arrays to store velocity increments at bottom and top boundaries; all functions of x
[delta_u_b,delta_u_t,delta_v_b,delta_v_t,delta_w_b,delta_w_t] = create_1d_arrays(nx)

# 1d arrays to store velocity increments at east and west boundaries; all functions of z
[delta_u_e,delta_u_w,delta_v_e,delta_v_w,delta_w_e,delta_w_w] = create_1d_arrays(nz)


#---------------------------------------------------------------------------
# compute initial conditions
#---------------------------------------------------------------------------
for k in range(nz):
	for i in range(nx):
		
		#---------------------------------------------------------------------------
		# map nested to global coordinates, evaluate exact soln at tn  for use as IC 
		#---------------------------------------------------------------------------
		t = 0.0
		X = x0 + x[i]
		Z = z0 + z[k]
		u[i,k] = parent_soln(X,Z,t,0)
		v[i,k] = parent_soln(X,Z,t,1)
		w[i,k] = parent_soln(X,Z,t,2)
		b[i,k] = parent_soln(X,Z,t,3)
					

#---------------------------------------------------------------------------
# time stepping loop would start here, time step counter starts at zero
#---------------------------------------------------------------------------
for n in range(0,nsteps):

	tn = n*dt 
	tnp1 = tn + dt
	
	if( np.mod(n,inc_write_netcdf)==0 ): 
		do_ncwrite = True
	else:
		do_ncwrite = False

	#---------------------------------------------------------------------------
	# evaluate the rhs arrays for u*, v*, w* and b  using [u,v,w,b] at tn
	#---------------------------------------------------------------------------
	u_x,u_z = grad(u,x,z,frac,LU_x,LU_z,Q,basis_functions)
	v_x,v_z = grad(v,x,z,frac,LU_x,LU_z,Q,basis_functions)
	w_x,w_z = grad(w,x,z,frac,LU_x,LU_z,Q,basis_functions)
	b_x,b_z = grad(b,x,z,frac,LU_x,LU_z,Q,basis_functions)

	if( nu > 0. ):
		g2u = divergence(u_x,u_z,x,z,frac,LU_x,LU_z,Q,basis_functions)
		g2v = divergence(v_x,v_z,x,z,frac,LU_x,LU_z,Q,basis_functions)
		g2w = divergence(w_x,w_z,x,z,frac,LU_x,LU_z,Q,basis_functions)
	if( kappa > 0. ):
		g2b = divergence(b_x,b_z,x,z,frac,LU_x,LU_z,Q,basis_functions)

	#  rhs include rotation, inertia, buoyancy and diffusion , no pressure gradients
	MM0 = idx[0]
	RHS[:,:,0,MM0] =  f*v  - (u*u_x + w*u_z) + nu*g2u       # rhs for u*
	RHS[:,:,1,MM0] = -f*u  - (u*v_x + w*v_z) + nu*g2v       # rhs for v*
	RHS[:,:,2,MM0] =    b  - (u*w_x + w*w_z) + nu*g2w       # rhs for w*
	RHS[:,:,3,MM0] = -N2*w - (u*b_x + w*b_z) + kappa*g2b    # rhs for b



	#---------------------------------------------------------------------------
	# integrate the rhs terms one time step to get [ustar,vstar,wstar,b(tnp1)]
	# cyclic array idx updated in time_step 
	#---------------------------------------------------------------------------
	U = [u,v,w,b]                                                   # fields at t_n
	[ustar,vstar,wstar,b], idx = time_step(n,dt,U,RHS,AB_order,idx)



	#------------------------------------------------------------------------------------
	# determine the corrections needed for normal components at the boundaries 
	#   simply insert the BVALS for b 
	#   simply insert the BVALS for v  (only because d/dy=0 and v not pressure projected)
	#------------------------------------------------------------------------------------
	#    EAST & WEST
	for k in range(nz):
		U = parent_soln(x0,z0+z[k],tnp1,0)
		delta_u_e[k] = U - ustar[0,k]
	
		U = parent_soln(x0+Lx,z0+z[k],tnp1,0)
		delta_u_w[k] = U - ustar[nx-1,k]
	
		b[0,k]    = parent_soln(x0,z0+z[k],tnp1,3)
		b[nx-1,k] = parent_soln(x0+Lx,z0+z[k],tnp1,3)
	
		v[0,k]    = parent_soln(x0,z0+z[k],tnp1,1)
		v[nx-1,k] = parent_soln(x0+Lx,z0+z[k],tnp1,1)

	#    BOTTOM & TOP	
	for i in range(nx):
		W = parent_soln(x0+x[i],z0,tnp1,2)
		delta_w_b[i] = W - wstar[i,0]
	
		W = parent_soln(x0+x[i],z0+Lz,tnp1,2)
		delta_w_t[i] = W - wstar[i,nz-1]
	
		b[i,0]    = parent_soln(x0+x[i],z0,tnp1,3)
		b[i,nz-1] = parent_soln(x0+x[i],z0+Lz,tnp1,3)
	
		v[i,0]    = parent_soln(x0+x[i],z0,tnp1,1)
		v[i,nz-1] = parent_soln(x0+x[i],z0+Lz,tnp1,1)


	#-------------------------------------------------------------------------------------
	# add these to normal components of ustar, supported by approximate step functions 
	#-------------------------------------------------------------------------------------		
	for k in range(nz):
		for i in range(nx):		
			ustar[i,k] = ustar[i,k] + delta_u_e[k]*step_x_e[i] + delta_u_w[k]*step_x_w[i] 		                        
			wstar[i,k] = wstar[i,k] + delta_w_b[i]*step_z_b[k] + delta_w_t[i]*step_z_t[k]
		                        		

	#------------------------------------------------------------------------------
	# determine the corrections needed for tangential components at the boundaries 
	#------------------------------------------------------------------------------
	for k in range(nz):

		# EAST BOUNDARY  n_hat in negative x direction
		V_x,W_x = parent_derivs(x0,z0+z[k],tnp1,'x')
	
		ii=0   # left end
		ddx_v = endpoint_deriv(vstar[:,k],x,ii) 
		ddn_v = -ddx_v
		delta_v_e[k] = gamma_x*( -V_x - ddn_v )
	
		ii=0
		ddx_w = endpoint_deriv(wstar[:,k],x,ii) 
		ddn_w = -ddx_w
		delta_w_e[k] = gamma_x*( -W_x - ddn_w )
	
		# WEST BOUNDARY  n_hat in positive x direction
		V_x,W_x = parent_derivs(x0+Lx,z0+z[k],tnp1,'x')
	
		ii=-1  # right end
		ddx_v = endpoint_deriv(vstar[:,k],x,ii)
		ddn_v = ddx_v
		delta_v_w[k] = gamma_x*( V_x - ddn_v )
	
		i=nx-1
		ddx_w = endpoint_deriv(wstar[:,k],x,ii)
		ddn_w = ddx_w
		delta_w_w[k] = gamma_x*( W_x - ddn_w )

	for i in range(nx):
	
		# BOTTOM BOUNDARY  n_hat in negative z direction
		U_z,V_z = parent_derivs(x0+x[i],z0,tnp1,'z')
	
		kk=0   # bottom
		ddz_u = endpoint_deriv(ustar[i,:],z,kk)
		ddn_u = -ddz_u
		delta_u_b[i] = gamma_z*( -U_z - ddn_u )
	
		ddz_v = endpoint_deriv(vstar[i,:],z,kk)
		ddn_v = -ddz_v
		delta_v_b[i] = gamma_z*( -V_z - ddn_v )
	
		# TOP BOUNDARY  n_hat in positive z direction
		U_z,V_z = parent_derivs(x0+x[i],z0+Lz,tnp1,'z')
	
		kk=-1   # top
		ddz_u = endpoint_deriv(ustar[i,:],z,kk)
		ddn_u = ddz_u
		delta_u_b[i] = gamma_z*( U_z - ddn_u )
	
		k=nz-1
		ddz_v = endpoint_deriv(vstar[i,:],z,kk)
		ddn_v = ddz_v
		delta_v_b[i] = gamma_z*( V_z - ddn_v )

	#-------------------------------------------------------------------------------------
	# add these to tangential components of ustar, supported by decaying exponentials 
	#  don't do this for v, but only because d/dy=0 and v is not pressure projected
	#-------------------------------------------------------------------------------------
	for k in range(nz):
		for i in range(nx):	
	
			ustar[i,k] = ustar[i,k] + delta_u_b[i]*exp_z_b[k] + delta_u_t[i]*exp_z_t[k]
		                        	
			#vstar[i,k] = vstar[i,k] + delta_v_e[k]*exp_x_e[i] + delta_v_w[k]*exp_x_w[i] \
			#                        + delta_v_b[i]*exp_z_b[k] + delta_v_t[i]*exp_z_t[k]
		                        
			wstar[i,k] = wstar[i,k] + delta_w_e[k]*exp_x_e[i] + delta_w_w[k]*exp_x_w[i] 		                        


	#-----------------------------------------------------------------------------
	# compute div u* = d/dx u* + d/dz w*  
	# assume no particular symmetries for u* and w* 
	#-----------------------------------------------------------------------------
	divstar = divergence(ustar,wstar,x,z,frac,LU_x,LU_z,Q,basis_functions)
	divmean = np.mean(divstar)
	#divstar = divstar - divmean     # doesn't seem to matter


	#---------------------------------------------------------------------------
	#  invert the div of u_s^* to generate an approximate solution for phi
	#  poisson equation has homogeneous Neumann conditions
	#---------------------------------------------------------------------------
	phi = poisson_invert_2d(divstar/dt,inversion_matrix)


	#---------------------------------------------------------------------------
	#  compute the gradient of phi
	#---------------------------------------------------------------------------
	phi_x,phi_z = grad(phi,x,z,frac,LU_x,LU_z,Q,basis_functions)


	#---------------------------------------------------------------------------
	#  do the pressure projection to get velocity components at tnp1
	#   (N.B. this is not a time integration via an Euler step)
	#---------------------------------------------------------------------------
	u = ustar - dt*phi_x
	w = wstar - dt*phi_z


	#---------------------------------------------------------------------------
	#  just to confirm that the projected velocity field is divergence free
	#---------------------------------------------------------------------------
	div = divergence(u,w,x,z,frac,LU_x,LU_z,Q,basis_functions)

		

	if( do_ncwrite ):
		#---------------------------------------------------------------------------
		# compute exact soln at tnp1
		#---------------------------------------------------------------------------
		u_exact=np.zeros_like(u); v_exact=np.zeros_like(v) 
		w_exact=np.zeros_like(w); b_exact=np.zeros_like(b)

		for k in range(nz):
			for i in range(nx):
		
				#---------------------------------------------------------------------------
				# map nested to global coordinates, evaluate exact soln at tn  for use as IC 
				#---------------------------------------------------------------------------
				X = x0 + x[i]
				Z = z0 + z[k]	
				u_exact[i,k] = parent_soln(X,Z,tnp1,0)
				v_exact[i,k] = parent_soln(X,Z,tnp1,1)
				w_exact[i,k] = parent_soln(X,Z,tnp1,2)
				b_exact[i,k] = parent_soln(X,Z,tnp1,3)
				
		U0 = A
		B0 = A*(k_iw/m_iw)*(N2/omega)
		rms_u = np.sqrt( np.sum( (u-u_exact)**2 )/(nx*nz) )/U0	
		rms_v = np.sqrt( np.sum( (v-v_exact)**2 )/(nx*nz) )/U0
		rms_w = np.sqrt( np.sum( (w-w_exact)**2 )/(nx*nz) )/U0
		rms_b = np.sqrt( np.sum( (b-b_exact)**2 )/(nx*nz) )/B0	
		print('RMS error in u, scaled by U0:  ',rms_u)
		print('RMS error in v, scaled by U0:  ',rms_v)
		print('RMS error in w, scaled by U0:  ',rms_w)
		print('RMS error in b, scaled by B0:  ',rms_b)

		#------------------------------------------------------------------------------
		# Output everything
		#------------------------------------------------------------------------------	
		fname = fname = "output/XZ_" + str(n).zfill(4) + ".nc"                                                      # construct the netcdf filename
		vars = [ustar,wstar,divstar,phi,u,v,w,b,u_exact,v_exact,w_exact,b_exact,div]                                # list of arrays with data
		varnames = ['ustar','wstar','div*','phi','u','v','w','b','u_exact','v_exact','w_exact','b_exact','div_u']   # list of variable names
		dims = ['m/s','m/s','1/s','m2/s2','m/s','m/s','m/s','m/s2','m/s','m/s','m/s','m/s2','1/s']                  # list of dimensions
		success = write_netcdf_2d(vars,x,z,tnp1,fname,varnames,dims)              # create and write the netcdf file	
