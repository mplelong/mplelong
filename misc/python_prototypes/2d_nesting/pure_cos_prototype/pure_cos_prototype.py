#!/usr/local/bin/python
"""
  Illustrate/explore the projection scheme for boundary forced problems
  limiting the numerics to straight cosine expansions with term by term
  differentiation.
    
  Numerics:
  (1) Use term by term differentiation of cos expansions of all variables
  (2) cosine-transform based fast solver for Poisson w/ homogeneous Neumann BCs
  (3) boundary information incorporated into the vector u* so that this is the
      elliptic equation we really need to solve
  (4) some additional "fixes" made to mitigate the changes we make to the problem
      by altering functional behavior near the boundaries and to compensate for
      inaccuracies in computing some quantities near the boundaries.
  
  --------------------------------------------------------
  Test problem -- exact solution in outer domain (X,Z):
  2D propagating internal wave mode: 
  --------------------------------------------------------
  	U = A cos(m*Z)*cos(k*X-omega*t) 
  	V = A*(f/omega)*cos(m*Z)*sin(k*X-omega*t) 
  	W = A*(k/m)*sin(m*Z)*sin(k*X-omega*t) 
  	B = -A*(k/m)(N^2/omega)*sin(m*Z)*cos(k*X-omega*t) 
  	with omega^2 = (k^2*N^2 + m^2*f^2)/(k^2 + m^2) 
  	m = (pi/H)*{1,2,3...}    k = (2pi/L)*{1,2,3...} 
  	NB  B=W=0 at z=0,H    flow is periodic in x with periodicity L
  
  --------------------------------------------------------	
  Governing equations in 2d nested domain (x,z):  
       w/ 1/rho_0 absorbed into definition of p
       b = -g/rho_0 * rho'   rho = rho_bar(z) + rho'
       N^2 = -(g/rho_0)*d/dz(rho_bar) = constant > 0
  -------------------------------------------------------- 
    u_x + w_z = 0     
  	u_t =  fv - (u u_x + w u_z) - p_x + nu Diff(u)           
  	v_t = -fu - (u v_x + w v_z) - p_y + nu Diff(v)   (w/ p_y = 0)
  	w_t =   b - (u w_x + w w_z) - p_z + nu Diff(w)             
  	b_t = - N^2 w - (u b_x + w b_z) + kappa Diff(b)           
    	  
  BCs:
   [u,v,w,b](t) and normal derivatives as necessary:
    specified by evaluating [U,V,W,B] on boundary of nested domain
  --------------------------------------------------------
  
  
  --------------------------------------------------------------------  
  Nested domain: 
  	x in [0,Lx]  z in [0,Lz]  offset from outer domain by x0 and z0
  	==> X = x0 + x  and  Z = z0 + z
  --------------------------------------------------------------------
  
  	
  2D array storage index convention:
  u[x,z] etc, number of gridpoints on closed interior domain = [nx,nz]
   
  variable ids = [0,1,2,3,4,5,6] for [u,v,w,b,eta,zeta,p] respectively
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from solver_utilities import create_inversion_matrix, poisson_invert_2d
from solver_utilities import create_diffusion_matrix, diffuse_2d
from solver_utilities import create_2d_variables, create_1d_arrays
from solver_utilities import apply_bcs, zero_derivs_near_boundaries, time_step,   \
                             add_nudging, divergence, div_ustar, grad, step, parent_soln, \
                             fill_boundary_vals, extrap_near_boundaries
from solver_utilities import add_body_force

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
nu = 1.e-3                # [m2/s]        viscosity
kappa = nu                # [m2/s]        diffusivity

k_iw = 2.*pi/L            # [1/m]         horizontal wavenumber of internal wave
m_iw = 1.*pi/H            # [1/m]         vertical wavenumber of internal wave

omega2 = (k_iw**2 * N2 + m_iw**2 * f2)/(k_iw**2 + m_iw**2)

omega = np.sqrt(omega2)   # [1/s]         internal wave frequency
IWPER = (2.*pi/omega)     # [s]
A = 0.01                  # [m/s]         amplitude of U
phase = pi/4.             # [1]           initial phase shift of wave mode
WAVE_PARAMS = [A,f,N2,omega,k_iw,m_iw,phase]
	

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#  parameters to define and control the nested simulation
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
Lx = L/5.                        # [m]      size of nested domain in x  10 km
Lz = H/5.                        # [m]      size of nested domain in z
x0 = 0.50*L				         # [m]      horizontal offset for origin of nested domain
z0 = 0.60*H                      # [m]      vertical offset for origin of nested domain

# x nested grid 	
nx = 513                         # [1] number of discrete grid points, basis functions
dx = Lx/(nx-1.)                  # [m] grid spacing
x = np.linspace(0,Lx,nx)         # [m] discrete grid points in [0,Lx]

# z nested grid 	
nz = 257                         # [1] number of discrete grid points, basis functions
dz = Lz/(nz-1.)                  # [m] grid spacing
z = np.linspace(0,Lz,nz)         # [m] discrete grid points in [0,Lz]

#------------------------------------------------------------
# parameters controlling the discrete time integration
#------------------------------------------------------------
nt = 1024                        # [1] time steps per wave period  unstable at 129, 257 seems ok
dt = IWPER/nt                    # [s] time step 
nsteps = nt+1
inc_write_netcdf = (nsteps-1)/64  
iplot = nsteps-1
AB_order = 4                     # 1 is flag for Euler. AB_order <= 4.

T_diff = 3*dt                    # diffusion time at Nyquist scale
p = 4                            # 1/2 order of diffusion operator
k_nyq = (np.pi/Lx)*(nx-1)
nu_star_x = k_nyq**(2*p)/T_diff  # save and use this value for momentum and scalars
k_nyq = (np.pi/Lz)*(nz-1)
nu_star_z = k_nyq**(2*p)/T_diff  # save and use this value for momentum and scalars


#------------------------------------------------------------
# numerical differentiation parameters
#------------------------------------------------------------
frac = 0.075                  # Fourier wavenumber filtering parameter

print( "domain aspect Lx/Lz  ",Lx/Lz)
print( "dx  in m             ",dx)
print( "dz  in m             ",dz)
print( "dt  in seconds       ",dt)
print( "cfl x  U dt/dx       ",A*dt/dx )
print( "cfl z  W dt/dz       ",A*(k_iw/m_iw)*dt/dz )
print( "diff cfl nu dt/dx^2  ",nu*dt/dx**2 )
print( "diff cfl nu dt/dz^2  ",nu*dt/dz**2 )
print( "IWPER  in hours      ",IWPER/3600.) 
print( "IPER   in hours:     ",IPER)
print( "BVPER  in mins:      ",BVPER*60)
print( "omega/f              ",omega/f)
print( "N/f                  ",np.sqrt(N2)/f)



#-------------------------------------------------------------------------------------
#     PRELIMINARY TASKS
#-------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
# set width of adjustment zones in grid points
# construct the approximate step functions I'll use near boundaries
#--------------------------------------------------------------------------
npts=2    # IW test problem works for 1 period w/ npts=2   
gamma = npts*dx
step_e,step_w = step(x,gamma)          # east, west
#--------------------------------------------------------------------------
gamma = 3*npts*dz
step_b,step_t = step(z,gamma)          # bottom, top
STEP = [step_e,step_w,step_b,step_t]   # pack up the step functions for ease of use
	
#----------------------------------------------------------------------------------
# create the (extended) inversion matrix  -1/(k^2+m^2)  for Fourier pressure solve 
# (will be applied in 2d wavenumber space to even extension of div ustar)
#----------------------------------------------------------------------------------
inversion_matrix = create_inversion_matrix(nx,nz,Lx,Lz,frac)

#--------------------------------------------------------------------------
# create the diffusion matrices for analytic integration 
# of diffusion terms in wavenumber space
#--------------------------------------------------------------------------
KAPPA=[nu_star_x,nu_star_z]   # different coeff for velocity in x and z directions
P = [p,p]                     # same half order of operator in x and z directions
diffusion_matrix_v = create_diffusion_matrix(nx,nz,Lx,Lz,frac,KAPPA,P,dt)

KAPPA=[nu_star_x,nu_star_z]   # different coeff for scalar b in x and z directions
P = [p,p]                     # same half order of operator in x and z directions
diffusion_matrix_b = create_diffusion_matrix(nx,nz,Lx,Lz,frac,KAPPA,P,dt)



#---------------------------------------------------------------------------------------
#  allocate arrays and fill with zeros
#---------------------------------------------------------------------------------------

# initialize some 2d arrays   (routine initializes 8 at a time, column major storage)
[u,v,w,b,ustar,vstar,wstar,divustar] = create_2d_variables(nx,nz)
phi=u.copy(); phi_x=u.copy(); phi_y=u.copy(); phi_z=u.copy(); div=u.copy() 

# initialize AB_order time levels of the 2d rhs for each of the 4 variables u, v, w and b
method = 'euler'  # initialize with 'euler', time stepping routine will update to match AB_order
nvars = 4         # u,v,w,b get time stepped varid=[0,1,2,3]
RHS = np.zeros((nx,nz,nvars,AB_order),float, order='F')     

#initialize the cyclic array for keeping track of previous time levels          
idx = [0,1,2,3,method]   # i.e. [MM0,MM1,MM2,MM3,method], will be toggled in the time_step routine

# some more arrays we'll need...
u_exact=np.zeros_like(u); v_exact=np.zeros_like(v) 
w_exact=np.zeros_like(w); b_exact=np.zeros_like(b); p_exact=np.zeros_like(b)

# 1d arrays to store boundary values at bottom and top boundaries; all 1d functions of x
[U_bot,V_bot,W_bot,B_bot,U_top,V_top,W_top,B_top] = create_1d_arrays(nx)
[U_z_bot,V_z_bot,W_z_bot,B_z_bot,U_z_top,V_z_top,W_z_top,B_z_top] = create_1d_arrays(nx)

# 1d arrays to store boundary values at east and west boundaries; all 1d functions of z
[U_east,V_east,W_east,B_east,U_west,V_west,W_west,B_west] = create_1d_arrays(nz)
[U_x_east,V_x_east,W_x_east,B_x_east,U_x_west,V_x_west,W_x_west,B_x_west] = create_1d_arrays(nz)

# pack things up for easy passing to functions
east_vals = [U_east,V_east,W_east,B_east]
west_vals = [U_west,V_west,W_west,B_west]
bot_vals =  [U_bot,V_bot,W_bot,B_bot]
top_vals =  [U_top,V_top,W_top,B_top]
BVALS = [east_vals,west_vals,bot_vals,top_vals]

east_derivs = [U_x_east,V_x_east,W_x_east,B_x_east] 
west_derivs = [U_x_west,V_x_west,W_x_west,B_x_west]
bot_derivs = [U_z_bot,V_z_bot,W_z_bot,B_z_bot]
top_derivs = [U_z_top,V_z_top,W_z_top,B_z_top] 
BDERIVS = [east_derivs,west_derivs,bot_derivs,top_derivs]


# location for plotting 1d slices
K = (nz-1)/4  ;  I = (nx-1)/4


#--------------------------------------------------------------------------------
# compute the initial conditions as a snapshot of the exact global  soln at t=0
#--------------------------------------------------------------------------------
for k in range(nz):
	for i in range(nx):		
		#----------------------------------------------------------------
		# map nested to global coordinates, evaluate exact soln at t=0 
		#----------------------------------------------------------------
		t = 0.0
		X = x0 + x[i]
		Z = z0 + z[k]
		u[i,k] = parent_soln(X,Z,t,0,WAVE_PARAMS)    # 4th arg is variable id
		v[i,k] = parent_soln(X,Z,t,1,WAVE_PARAMS)
		w[i,k] = parent_soln(X,Z,t,2,WAVE_PARAMS)
		b[i,k] = parent_soln(X,Z,t,3,WAVE_PARAMS)



#---------------------------------------------------------------------------
# time stepping loop starts here, time step counter starts at zero
#---------------------------------------------------------------------------
for n in range(0,nsteps):

	tn = n*dt 
	tnp1 = tn + dt
	
	if( np.mod(n,inc_write_netcdf)==0 ): 
		do_ncwrite = True
	else:
		do_ncwrite = False
	
	#--------------------------------------------------------------------------
	# Adjust the starting point of each time step so that normal derivatives of 
	# all fields vanish at the boundaries. Derivatives are then computed using
	# term by term differentiation of cos series expansions of all dependent  
	# variables and intermediate quantities.
	#---------------------------------------------------------------------------
	dir='xz'
	u = zero_derivs_near_boundaries(u,x,z,npts,dir)
	v = zero_derivs_near_boundaries(v,x,z,npts,dir)
	w = zero_derivs_near_boundaries(w,x,z,npts,dir)
	b = zero_derivs_near_boundaries(b,x,z,npts,dir)					

	
	if( do_ncwrite or n==iplot ):	
		#--------------------------------------------------------------
		# compute exact soln at tnp1, i.e. at end of this time step
		#--------------------------------------------------------------
		for k in range(nz):
			for i in range(nx):		
				X = x0 + x[i]
				Z = z0 + z[k]	
				u_exact[i,k] = parent_soln(X,Z,tnp1,0,WAVE_PARAMS)
				v_exact[i,k] = parent_soln(X,Z,tnp1,1,WAVE_PARAMS)
				w_exact[i,k] = parent_soln(X,Z,tnp1,2,WAVE_PARAMS)
				b_exact[i,k] = parent_soln(X,Z,tnp1,3,WAVE_PARAMS)
				p_exact[i,k] = parent_soln(X,Z,tnp1,6,WAVE_PARAMS)		

	#---------------------------------------------------------------------------
	# evaluate the rhs arrays for u*, v*, w* and b  using [u,v,w,b] at tn
	#---------------------------------------------------------------------------
	nvals=3*npts
	u_x,u_z = grad(u,x,z,frac,nvals)
	v_x,v_z = grad(v,x,z,frac,nvals)
	w_x,w_z = grad(w,x,z,frac,nvals)
	b_x,b_z = grad(b,x,z,frac,nvals)
	
	#-------------------------------------------------------------------------------
	#  rhs include rotation, inertia, buoyancy; no pressure gradients, no diffusion
	#-------------------------------------------------------------------------------
	MM0 = idx[0]
	RHS[:,:,0,MM0] =  f*v   - (u*u_x + w*u_z)        # rhs for u*
	RHS[:,:,1,MM0] = -f*u   - (u*v_x + w*v_z)        # rhs for v*
	RHS[:,:,2,MM0] =  b     - (u*w_x + w*w_z)        # rhs for w*   
	RHS[:,:,3,MM0] = -N2*w  - (u*b_x + w*b_z)        # rhs for b	
	
	#------------------------------------------------------------------------------------
	#  Get the boundary values at t=tnp1 and store in BVALS and BDERIVS. 
	#  Generally these are externally supplied vals. 
	#  Here we just evaluate and store the t=tnp1 outer solns. 
	#------------------------------------------------------------------------------------
	BVALS,BDERIVS = fill_boundary_vals(x,z,tnp1,x0,z0,WAVE_PARAMS,BVALS,BDERIVS)
	
	
	#---------------------------------------------------------------------------
	# Add the nudging terms relaxing toward externally supplied boundary values
	# I can't get this to help for this problem in any way.
	#---------------------------------------------------------------------------
	#U=[u,v,w,b]
	#RHS = add_nudging(RHS,U,x,z,BVALS,STEP,dt,MM0) 


	#------------------------------------------------------------------------------------
	# Integrate the explicit rhs terms one step to get [ustar,vstar,wstar,b(tnp1)].
	# The cyclic array idx is updated in time_step. 
	#------------------------------------------------------------------------------------
	U = [u,v,w,b]   # fields at t=t_n
	[ustar,vstar,wstar,b],idx = time_step(n,dt,U,RHS,AB_order,idx)
			
		
	if( n==iplot ):
		plt.plot( x/Lx,   ustar[:,K],'k' )
		plt.plot( x/Lx, u_exact[:,K],'r' )
		plt.plot( x/Lx,   wstar[:,K],'k' )
		plt.plot( x/Lx, w_exact[:,K],'b' )
		plt.xlabel('x/Lx')
		plt.ylabel('')
		plt.savefig('ustar_hat_of_x.eps',dpi=300,bb_inches='tight')
		plt.close()
				
		plt.plot( wstar[I,:],z/Lz,'k' )
		plt.plot( w_exact[I,:],z/Lz,'b' )
		plt.plot( ustar[I,:],z/Lz,'k' )
		plt.plot( u_exact[I,:],z/Lz,'r' )
		plt.ylabel('z/Lz')
		plt.xlabel('')
		plt.savefig('ustar_hat_of_z.eps',dpi=300,bb_inches='tight')
		plt.close()
	
	
	#-----------------------------------------------------------------------------
	# Compute div {\hat u*} = d/dx {\hat u*} + d/dz {\hat w*}.  
	#  The divergence is augmented with normal velocities at the boundaries
	#  so that the correct elliptic problem with homogeneous Neumann BCs
	#  for pressure is set up.
	#-----------------------------------------------------------------------------
	[divstar,ustar,wstar] = div_ustar(ustar,wstar,x,z,frac,npts,BVALS,STEP)	
	
			
	if( n==iplot ):
		plt.plot( x/Lx, (Lx/A)* divstar[:,K],'k' )
		plt.xlabel('x/Lx')
		plt.ylabel('div u* x (Lx/A)')
		plt.savefig('divustarhat_of_x.eps',dpi=300,bb_inches='tight')
		plt.close()
		
		plt.plot( (Lx/A)* divstar[I,:],z/Lz,'k' )
		plt.ylabel('z/Lz')
		plt.xlabel('div u* x (Lx/A)')
		plt.savefig('divustarhat_of_z.eps',dpi=300,bb_inches='tight')
		plt.close()

	#---------------------------------------------------------------------------
	#  Poisson-invert the div of {\hat u*} to obtain the pressure variable phi
	#  that satisfies the Poisson equation with homogeneous Neumann conditions.
	#  The solution is wavenumber filtered using the parameter frac.
	#---------------------------------------------------------------------------
	phi = poisson_invert_2d(divstar/dt,inversion_matrix)
	
	if( n==iplot ):
		plt.plot( x/Lx, phi[:,K],'k' )
		plt.plot( x/Lx, p_exact[:,K],'r' )
		plt.xlabel('x/Lx')
		plt.ylabel('P_hat, p_n+1')
		plt.savefig('P_hat_of_x.eps',dpi=300,bb_inches='tight')
		plt.close()
		
		plt.plot( phi[I,:], z/Lz,'k' )
		plt.plot( p_exact[I,:],z/Lz,'r' )
		plt.ylabel('z/Lz')
		plt.xlabel('P_hat, p_n+1')
		plt.savefig('P_hat_of_z.eps',dpi=300,bb_inches='tight')
		plt.close()

	#---------------------------------------------------------------------------
	#  Compute the gradient of phi.
	#---------------------------------------------------------------------------
	nvals=0  # no fd replacements near edges
	phi_x,phi_z = grad(phi,x,z,frac,nvals)

	if( n==iplot ):
		B0 = A*(k_iw/m_iw)*(N2/omega)
		plt.plot( x/Lx, phi_x[:,K]/B0,'r',x/Lx, phi_z[:,K]/B0,'b' )
		plt.xlabel('x/Lx')
		plt.ylabel('')
		plt.grid()
		plt.savefig('grad_P_hat_of_x.eps',dpi=300,bb_inches='tight')
		plt.close()
		
		plt.plot( phi_x[I,:]/B0, z/Lz,'r',phi_z[I,:]/B0, z/Lz,'b' )
		plt.ylabel('z/Lz')
		plt.xlabel('')
		plt.grid()
		plt.savefig('grad_P_hat_of_z.eps',dpi=300,bb_inches='tight')
		plt.close()

	#---------------------------------------------------------------------------
	#  Do the pressure projection to get velocity components at t=tnp1.
	#  (N.B. this is not a time integration via an Euler step)
	#---------------------------------------------------------------------------
	u = ustar - dt*phi_x
	v = vstar
	w = wstar - dt*phi_z
	
	#---------------------------------------------------------------------------
	#  Apply some fixes to clean up behavior near the boundaries.
	#  apply_bcs              uses externally supplied normal derivs at boundaries 
	#                         and extrapolates inward across BL with a small floating offset
	#
	#  extrap_near_boundaries uses the computed solution "just inside" the BL
	#                         and extrapolates out toward the boundary
	#
	#  my estimate of the BL width is ~3*npts = 6 points in this case
	#
	#  In one or a few time step tests, both schemes work well but over many steps
	#  boundary instabilities build up. Alternating between the schemes prevents
	#  this. 
	#---------------------------------------------------------------------------
	if( np.mod(n,2)==0 ):
		[u,v,w,b] = apply_bcs(u,v,w,b,x,z,BVALS,BDERIVS,npts)    # info to interior
	else:
		dir='xz'
		xx = 1.5    # works for 1 per w xx=2,1.5
		u = extrap_near_boundaries(u,x,z,xx*npts,dir)            # info from interior
		v = extrap_near_boundaries(v,x,z,xx*npts,dir)
		w = extrap_near_boundaries(w,x,z,xx*npts,dir)
		b = extrap_near_boundaries(b,x,z,xx*npts,dir)
	
	#---------------------------------------------------------------------------
	#  exact integration of diffusion
	#---------------------------------------------------------------------------
	u = diffuse_2d(u,diffusion_matrix_v)
	v = diffuse_2d(v,diffusion_matrix_v)
	w = diffuse_2d(w,diffusion_matrix_v)
	b = diffuse_2d(b,diffusion_matrix_b)
	
	
	
	#--------------------------------------------------------------------------------
	#  The computed solutions are our attempt at the true solutions that don't have
	#  artificially altered zero-normal-derivative behavior at the boundaries.
	#  N.B.:  to compute these solutions, we tweaked the starting fields near
	#         the boundaries, i.e. the time step algorithm relies on these initial
	#         tweaks ==> this becomes the first step of each time step. 
	#--------------------------------------------------------------------------------
	
					

	if( do_ncwrite ):	
		#---------------------------------------------------------------------------
		#  Just to confirm that the projected velocity field is divergence free.
		#---------------------------------------------------------------------------
		nvals=3*npts
		div = divergence(u,w,x,z,frac,nvals)
					
		U0 = A
		W0 = A*(k_iw/m_iw)
		B0 = A*(k_iw/m_iw)*(N2/omega)
		rms_u = np.sqrt( np.sum( (u-u_exact)**2 )/(nx*nz) )/U0	
		rms_v = np.sqrt( np.sum( (v-v_exact)**2 )/(nx*nz) )/U0
		rms_w = np.sqrt( np.sum( (w-w_exact)**2 )/(nx*nz) )/W0
		rms_b = np.sqrt( np.sum( (b-b_exact)**2 )/(nx*nz) )/B0	
		print('RMS error in u, scaled by U0:  ',rms_u)
		print('RMS error in v, scaled by U0:  ',rms_v)
		print('RMS error in w, scaled by W0:  ',rms_w)
		print('RMS error in b, scaled by B0:  ',rms_b)

		#------------------------------------------------------------------------------
		# Output everything
		#------------------------------------------------------------------------------	
		fname = fname = "output/XZ_" + str(n).zfill(4) + ".nc"                                                      # construct the netcdf filename
		vars = [ustar,wstar,divstar,phi,u,v,w,b,phi_x,phi_z,u_exact,v_exact,w_exact,b_exact,div]                                # list of arrays with data
		varnames = ['ustar','wstar','div*','phi','u','v','w','b','phi_x','phi_z','u_exact','v_exact','w_exact','b_exact','div_u']   # list of variable names
		dims = ['m/s','m/s','1/s','m2/s2','m/s','m/s','m/s','m/s2','m/s2','m/s2','m/s','m/s','m/s','m/s2','1/s']                  # list of dimensions
		success = write_netcdf_2d(vars,x,z,tnp1,fname,varnames,dims)                                                # create and write the netcdf file	


