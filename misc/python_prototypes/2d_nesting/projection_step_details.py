#!/usr/local/bin/python
"""
  Illustrate/explore the projection scheme for boundary forced problems
  using the nudging approach
    
  Numerics:
  (1) Use Bernoulli-Cos method for derivatives, no assumed symmetries 
  (2) cosine-transform based fast solver for Poisson w/ homogeneous Neumann BCs
  (3) u* vector updated with boundary-normal "jump" conditions
  
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
  	u_t =  fv - (u u_x + w u_z) + nu Laplacian(u) - p_x          
  	v_t = -fu - (u v_x + w v_z) + nu Laplacian(v)
  	w_t =   b - (u w_x + w w_z) + nu Laplacian(w) - p_z            
  	b_t = - N^2 w - (u b_x + w b_z) + kappa Laplacian(b)           
  	eta_t  = u
  	zeta_t = w
    	  
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
from solver_utilities import prepare_BernoulliCosine, evaluate_basis_functions
from solver_utilities import create_2d_variables, create_1d_arrays
from solver_utilities import time_step, divergence, grad, grad2
from solver_utilities import construct_psi_u, construct_psi_w, apply_bcs, add_body_force

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


def parent_soln(X,Z,t,id):
	# function to evaluate parent soln at X,Z,t  ; 
	# here an IW mode with wave parameters available globally
	# id = [0,1,2,3,6] for [U,V,W,B,P] respectively
	argz = m_iw*Z
	argx = k_iw*X - omega*t + phase
	if(id==0): ans = A*np.cos(argz)*np.cos(argx) # U
	if(id==1): ans = A*(f/omega)*np.cos(argz)*np.sin(argx) # V    
	if(id==2): ans = A*(k_iw/m_iw)*np.sin(argz)*np.sin(argx) # W 
	if(id==3): ans = -A*(k_iw/m_iw)*(N2/omega)*np.sin(argz)*np.cos(argx) # B
	if(id==6): ans = -A*(1./(k_iw*omega))*(f2-omega2)*np.cos(argz)*np.cos(argx) # P	
	return ans

def parent_derivs(X,Z,t,dir):
	# function to evaluate certain derivs of parent soln at X,Z,t  ; 
	# here an IW mode with wave parameters available globally
	# dir = 'x','y' or 'z'
	# used for normal derivs of tangential vels at boundaries
	argz = m_iw*Z
	argx = k_iw*X - omega*t + phase	
	if(dir=='x'):    # return W_x,W_xx  (in this order)
		ans1 = k_iw*A*(k_iw/m_iw)*np.sin(argz)*np.cos(argx) # W_x
		ans2 = -k_iw*k_iw*A*(k_iw/m_iw)*np.sin(argz)*np.sin(argx) # W_xx
		ans3 = 0.
		ans4 = 0.		
	if(dir=='y'):    # d/dy = 0 here 
		ans1 = 0.
		ans2 = 0.
		ans3 = 0.
		ans4 = 0
	if(dir=='z'):    # return U_z,U_zz,W_z,W_zz
		ans1 = -m_iw*A*np.sin(argz)*np.cos(argx)                   # U_z
		ans2 = -m_iw*m_iw*A*np.cos(argz)*np.cos(argx)              # U_zz
		ans3 = m_iw*A*(k_iw/m_iw)*np.cos(argz)*np.sin(argx)        # W_z
		ans4 = -m_iw*m_iw*A*(k_iw/m_iw)*np.sin(argz)*np.sin(argx)  # W_zz
	return ans1,ans2,ans3,ans4

	

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
nsteps = 1 # nt+1
inc_write_netcdf = 1 # (nsteps-1)/64  
iplot = nsteps-1
AB_order = 4                     # 1 is flag for Euler. AB_order <= 4.

#------------------------------------------------------------
# numerical differentiation parameters
#------------------------------------------------------------
Q = 5
frac = 0.05  # Fourier wavenumber filtering parameter

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
# create the (extended) inversion matrix  -1/k2  for Fourier pressure solve 
# (actually applied in 2d wavenumber space to even extended div ustar)
#--------------------------------------------------------------------------
inversion_matrix = create_inversion_matrix(nx,nz,Lx,Lz,frac)

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
#  initialize arrays
#---------------------------------------------------------------------------------------

# initialize some 2d arrays with zeros
[u,v,w,b,eta,zeta,ustar,vstar,wstar,g2u,g2v,g2w,g2b] = create_2d_variables(nx,nz)

# initialize 4 time levels of the 2d rhs for each of the 4 variables u, v, w and b
method = 'euler'  # always initialize with 'euler', time stepping routine will increment to match AB_order
nvars = 6         # u,v,w,b,eta and zeta get time stepped varid=[0,1,2,3]
RHS = np.zeros((nx,nz,nvars,AB_order),float)     

#initialize the cyclic array for defining previous time levels          
idx = [0,1,2,3,method]   # i.e. [MM0,MM1,MM2,MM3,method], will be toggled in the time_step routine

# 1d arrays to store velocity increments at bottom and top boundaries; all functions of x
[U_bot,V_bot,W_bot,B_bot,U_z_bot,U_zz_bot,W_z_bot,W_zz_bot] = create_1d_arrays(nx)
[U_top,V_top,W_top,B_top,U_z_top,U_zz_top,W_z_top,W_zz_top] = create_1d_arrays(nx)


# 1d arrays to store boundary values at east and west boundaries; all functions of z
[U_east,V_east,W_east,B_east,W_x_east,W_xx_east,U_x_east,U_xx_east] = create_1d_arrays(nz)
[U_west,V_west,W_west,B_west,W_x_west,W_xx_west,U_x_west,U_xx_west] = create_1d_arrays(nz)

u_exact=np.zeros_like(u); v_exact=np.zeros_like(v) 
w_exact=np.zeros_like(w); b_exact=np.zeros_like(b); p_exact=np.zeros_like(b)


# locations for plotting slices
K = int((nz-1)/4)  ;  I = int((nx-1)/4)


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
		u[i,k] = parent_soln(X,Z,t,0)    # last arg is variable id
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
	
	
	if( do_ncwrite or n==iplot ):	
		#--------------------------------------------------------------
		# compute exact soln at tnp1, i.e. at end of this time step
		#--------------------------------------------------------------
		for k in range(nz):
			for i in range(nx):		
				X = x0 + x[i]
				Z = z0 + z[k]	
				u_exact[i,k] = parent_soln(X,Z,tnp1,0)
				v_exact[i,k] = parent_soln(X,Z,tnp1,1)
				w_exact[i,k] = parent_soln(X,Z,tnp1,2)
				b_exact[i,k] = parent_soln(X,Z,tnp1,3)
				p_exact[i,k] = parent_soln(X,Z,tnp1,6)

	
	#------------------------------------------------------------------------------------
	#  get the required boundary information at tnp1 and store in 1d arrays 
	#  generally these are externally supplied vals, here; evaluate the tnp1 solns 
	#------------------------------------------------------------------------------------	
	for k in range(nz):	     #    save vals at EAST & WEST boundaries
		Z = z0+z[k]		
		# get the velocity components
		id = 0           
		U_east[k] = parent_soln(x0   ,Z,tnp1,id)
		U_west[k] = parent_soln(x0+Lx,Z,tnp1,id)
		id = 1           
		V_east[k] = parent_soln(x0   ,Z,tnp1,id)
		V_west[k] = parent_soln(x0+Lx,Z,tnp1,id)
		id = 2           
		W_east[k] = parent_soln(x0   ,Z,tnp1,id)
		W_west[k] = parent_soln(x0+Lx,Z,tnp1,id)
		# get the buoyancy
		id = 3           
		B_east[k] = parent_soln(x0   ,Z,tnp1,id)
		B_west[k] = parent_soln(x0+Lx,Z,tnp1,id)		
		# get the normal derivatives of the tangential velocities
		dir = 'x'
		W_x_east[k],W_xx_east[k],xx,yy = parent_derivs(x0   ,Z,tnp1,dir)
		W_x_west[k],W_xx_west[k],xx,yy = parent_derivs(x0+Lx,Z,tnp1,dir)
		
		
	for i in range(nx):	    #    save vals at BOTTOM & TOP boundaries
		X = x0 + x[i]		
		# get the velocity components
		id = 0
		U_bot[i] = parent_soln(X,z0   ,tnp1,id)
		U_top[i] = parent_soln(X,z0+Lz,tnp1,id)
		id = 1
		V_bot[i] = parent_soln(X,z0   ,tnp1,id)
		V_top[i] = parent_soln(X,z0+Lz,tnp1,id) 
		id = 2
		W_bot[i] = parent_soln(X,z0   ,tnp1,id)
		W_top[i] = parent_soln(X,z0+Lz,tnp1,id)
		# get the buoyancy
		id = 3           
		B_bot[i] = parent_soln(X,z0   ,tnp1,id)
		B_top[i] = parent_soln(X,z0+Lz,tnp1,id)		
		# get the normal derivatives of the tangential velocities
		dir = 'z'
		U_z_bot[i],U_zz_bot[i],W_z_bot[i],W_zz_bot[i] = parent_derivs(X,z0   ,tnp1,dir)
		U_z_top[i],U_zz_top[i],W_z_top[i],W_zz_top[i] = parent_derivs(X,z0+Lz,tnp1,dir)
		
		

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

	#  rhs include rotation, inertia, buoyancy and diffusion; no pressure gradients
	MM0 = idx[0]
	RHS[:,:,0,MM0] =  f*v   - (u*u_x + w*u_z) + nu*g2u       # rhs for u*
	RHS[:,:,1,MM0] = -f*u   - (u*v_x + w*v_z) + nu*g2v       # rhs for v*
	RHS[:,:,2,MM0] =  b     - (u*w_x + w*w_z) + nu*g2w       # rhs for w*   
	RHS[:,:,3,MM0] = -N2*w  - (u*b_x + w*b_z) + kappa*g2b    # rhs for b	
	RHS[:,:,4,MM0] =  u                                      # d/dt eta  = u
	RHS[:,:,5,MM0] =  w                                      # d/dt zeta = w

	#  add an oscillating "point source" to the eqn for w
	#Ampl = 1.e-4    # m/s2
	#sigma = .3*np.sqrt(N2)
	#RHS[:,:,2,MM0] = add_body_force(sigma,tn,x,z,Ampl,RHS[:,:,2,MM0])

	#---------------------------------------------------------------------------
	# integrate the rhs terms one time step to get [ustar,vstar,wstar,b(tnp1)]
	# cyclic array idx updated in time_step 
	#---------------------------------------------------------------------------
	U = [u,v,w,b,eta,zeta]   # fields at t=t_n
	[ustar,vstar,wstar,b,eta,zeta],idx = time_step(n,dt,U,RHS,AB_order,idx)
			

	#---------------------------------------------------------------------------
	# construct the function psi_u containing the boundary information on u 
	#---------------------------------------------------------------------------
	BVALS=[U_east,U_west,U_z_bot,U_z_top,ustar,wstar]
	psi_u = construct_psi_u(x,z,BVALS)	
	
	#---------------------------------------------------------------------------
	# construct the function psi_w containing the boundary information on w 
	#---------------------------------------------------------------------------
	BVALS=[W_bot,W_top,W_x_east,W_x_west,ustar,wstar]
	psi_w = construct_psi_w(x,z,BVALS)
	
	#------------------------------------------------------------------------------
	# update u* and w* to bring boundary inhomogeneity into the interior 
	# (in my notes, these modified variables are called u_star_hat and w_star_hat)
	#------------------------------------------------------------------------------
	ustar = ustar + psi_u
	wstar = wstar + psi_w
	
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
	# compute div {\hat u*} = d/dx {\hat u*} + d/dz {\hat w*}  
	# assume no particular symmetries for {\hat u*} and {\hat w*} 
	#-----------------------------------------------------------------------------
	divstar = divergence(ustar,wstar,x,z,frac,LU_x,LU_z,Q,basis_functions)
		
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
	#  Invert the div of {\hat u*} to obtain the pressure variable phi.
	#  phi satisfies the Poisson equation with homogeneous Neumann conditions.
	#  The solution is wavenumber filtered with parameter frac.
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
	#  compute the gradient of phi
	#---------------------------------------------------------------------------
	Q0 = 0              # use standard cos series differentiation
	phi_x,phi_z = grad(phi,x,z,frac,LU_x,LU_z,Q0,basis_functions)

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
	#  Do the pressure projection to get velocity components at tnp1.
	#  (N.B. this is not a time integration via an Euler step)
	#---------------------------------------------------------------------------
	u = ustar - dt*phi_x
	v = vstar
	w = wstar - dt*phi_z
	
	#---------------------------------------------------------------------------
	#  Impose the remaining BCS (scalar b, tangential velocity components)
	#---------------------------------------------------------------------------
	BVALS_EW = [U_east,U_west,V_east,V_west,W_east,W_west,B_east,B_west,W_x_east,W_x_west,W_xx_east,W_xx_west]
	BVALS_BT = [U_bot,U_top,V_bot,V_top,W_bot,W_top,B_bot,B_top,U_z_bot,U_z_top,U_zz_bot,U_zz_top,W_z_bot,W_z_top,W_zz_bot,W_zz_top]
	BVALS = [BVALS_EW,BVALS_BT]
	[u,v,w,b] = apply_bcs(u,v,w,b,x,z,BVALS)
				
	if( n==iplot ):
		iend=nx	
		inc = 4
		plt.plot(x/Lx,u_exact[:,K],'r')
		plt.plot(x[0:nx:inc]/Lx,u[0:nx:inc,K],'k.',markersize=4)
		plt.plot(x/Lx,w_exact[:,K],'b')
		plt.plot(x[0:nx:inc]/Lx,w[0:nx:inc,K],'k.',markersize=4)
		plt.xlabel(r"$x/L_x$")
		plt.savefig('compare_vs_x.eps',dpi=300,bb_inches='tight')
		plt.close()
	
		inc = 4
		plt.plot(w_exact[I,:],z/Lz,'b')
		plt.plot(w[I,0:nz:inc],z[0:nz:inc]/Lz,'k.',markersize=3)
		plt.plot(u_exact[I,:],z/Lz,'r')
		plt.plot(u[I,0:nz:inc],z[0:nz:inc]/Lz,'k.',markersize=3)
		plt.ylabel(r"$z/L_z$")
		plt.savefig('compare_vs_z.eps',dpi=300,bb_inches='tight')
		plt.close()
				

	if( do_ncwrite ):	
		#---------------------------------------------------------------------------
		#  Just to confirm that the projected velocity field is divergence free.
		#---------------------------------------------------------------------------
		div = divergence(u,w,x,z,frac,LU_x,LU_z,Q,basis_functions)
					
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


