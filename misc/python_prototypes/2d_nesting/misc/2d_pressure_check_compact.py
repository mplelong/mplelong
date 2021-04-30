#!/usr/local/bin/python
"""
  HERE WE TAKE EXACTLY ONE EULER TIME STEP (with tiny step size to approximate a smarter scheme)
  ===>  seems like good results in arbitrary interior subdomains
  
  
  Illustrate solving for 2d flow in a nested domain with open boundaries.
  BCs are obtained from analytic expressions for an internal wave mode propagating 
  through an arbitrary interior x-z subdomain of interest. The full domain has a 
  solid bottom and rigid lid w/ free-slip, no normal flow bcs on the outer soln.
  
  Numerics:
  (1) spectral-like compact schemes for derivatives, no symmetries assumed.
  (2) cosine-transform based fast solver for Poisson w/ homogeneous Neumann BCs
  (3) BCs on normal flow components embedded into source equation for pressure
  (4) BCs on normal dewrivatives of tangential velocity components also embedded
  
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
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from KW_FUNCTIONS import write_netcdf_2d,create_inversion_matrix,poisson_invert_2d,      \
                         build_compact_matrices_1,build_compact_matrices_2,compact_deriv,\
                         divergence,grad,grad2
#-----------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#  parameters for the environment and the outer solution which is used to
#  supply initial and time dependent boundary conditions
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
pi = np.pi
H = 4000.                 # [m]           full water depth
L = 1.e5                  # [m]           horizontal scale = 100 km
N = 2.0e-3                # [1/s]         buoyancy frequency  2
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
	if(id==1): ans = A*np.cos(argz)*np.cos(argx) # U
	if(id==2): ans = A*(f/omega)*np.cos(argz)*np.sin(argx) # V 
	if(id==3): ans = A*(k_iw/m_iw)*np.sin(argz)*np.sin(argx) # W 
	if(id==4): ans = -A*(k_iw/m_iw)*(N2/omega)*np.sin(argz)*np.cos(argx) # B
	if(id==5): ans = -A*(1./(k_iw*omega))*(f2-omega2)*np.cos(argz)*np.cos(argx) # P	
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
gamma_x = 8.*dx                     # [m] width of near boundary transition region


# z nested grid arrays	
nz = 129                            # [1] number of discrete grid points, basis functions
dz = Lz/(nz-1.)                     # [m] grid spacing
z = np.linspace(0,Lz,nz)            # [m] discrete grid points in [0,Lz]
gamma_z = 8.*dz                     # [m] width of near boundary transition region


# parameters controlling the discrete integration
nt = 1024                           # [1] time steps per wave period
dt = IWPER/nt                       # [s] time step 



#-------------------------------------------------------------------------------------
#     PRELIMINARY TASKS
#-------------------------------------------------------------------------------------


#--------------------------------------------------------------------------
# create the (extended) inversion matrix  -1/k2  for Fourier pressure solve 
# (actually applied in 2d wavenumber space to even extended div ustar)
#--------------------------------------------------------------------------
inversion_matrix = create_inversion_matrix(nx,nz,Lx,Lz)

#--------------------------------------------------------------------------
# create and factor the compact difference matrices needed for 1st and 2nd
# derivs in x and z directions
#--------------------------------------------------------------------------
lu_x_1,piv_x_1,B_x_1 = build_compact_matrices_1(nx,Lx)
lu_z_1,piv_z_1,B_z_1 = build_compact_matrices_1(nz,Lz)

lu_x_2,piv_x_2,B_x_2 = build_compact_matrices_2(nx,Lx)
lu_z_2,piv_z_2,B_z_2 = build_compact_matrices_2(nz,Lz)


#---------------------------------------------------------------------------------------
#  construct some useful arrays
#---------------------------------------------------------------------------------------
#  step_x is an array that is 1 at the x boundaries but vanishes in the interior
step_x_e = np.zeros_like(x); step_x_w = np.zeros_like(x)
#  exp_x is an array that decays from E/W boundaries with slope 1/gamma_x
exp_x_e = np.zeros_like(x); exp_x_w = np.zeros_like(x)
for i in range(nx):
	step_x_e[i] = np.exp(-((x[i]-0.)/gamma_x)**2) # approximate, discrete step function
	step_x_w[i] = np.exp(-((x[i]-Lx)/gamma_x)**2)
	exp_x_e[i] = np.exp(-((x[i]-0.)/gamma_x))
	exp_x_w[i] = np.exp(-((Lx-x[i])/gamma_x))
	
#  step_z is an array that is 1 at the z boundaries but vanishes in the interior
step_z_b = np.zeros_like(z); step_z_t = np.zeros_like(z)
#  exp_z is an array that decays from b/t boundaries with slope 1/gamma_z
exp_z_b = np.zeros_like(z); exp_z_t = np.zeros_like(z)
for k in range(nz):
	step_z_b[k] = np.exp(-((z[k]-0.)/gamma_z)**2)
	step_z_t[k] = np.exp(-((z[k]-Lz)/gamma_z)**2)
	exp_z_b[k] = np.exp(-((z[k]-0.)/gamma_z))
	exp_z_t[k] = np.exp(-((Lz-z[k])/gamma_z))

# initialize some arrays so they can be accessed elementwise
u = np.zeros((nx,nz),dtype=float) ; v = np.zeros((nx,nz),dtype=float)
w = np.zeros((nx,nz),dtype=float) ; b = np.zeros((nx,nz),dtype=float)
ustar = np.zeros((nx,nz),dtype=float) ; vstar = np.zeros((nx,nz),dtype=float)
wstar = np.zeros((nx,nz),dtype=float) ; 

# to hold corrections for normal components, i.e. for ustar at east and west boundaries
delta_u_e = np.zeros_like(z); delta_u_w = np.zeros_like(z);    # east and west
delta_w_b = np.zeros_like(x); delta_w_t = np.zeros_like(x);    # top and bottom

# to hold corrections for tangential components, i.e. for wstar at east and west boundaries
delta_w_e = np.zeros_like(z); delta_w_w = np.zeros_like(z);    # east and west
delta_v_e = np.zeros_like(z); delta_v_w = np.zeros_like(z);    # east and west
delta_u_b = np.zeros_like(x); delta_u_t = np.zeros_like(x);    # top and bottom
delta_v_b = np.zeros_like(x); delta_v_t = np.zeros_like(x);    # top and bottom










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
		u[i,k] = parent_soln(X,Z,t,1)
		v[i,k] = parent_soln(X,Z,t,2)
		w[i,k] = parent_soln(X,Z,t,3)
		b[i,k] = parent_soln(X,Z,t,4)
					




tn=0.0; 
tnp1=tn+dt
#---------------------------------------------------------------------------
# construct the rhs arrays for u,v,w and b
#---------------------------------------------------------------------------
u_x,u_z = grad(u,lu_x_1,piv_x_1,B_x_1,lu_z_1,piv_z_1,B_z_1)
v_x,v_z = grad(v,lu_x_1,piv_x_1,B_x_1,lu_z_1,piv_z_1,B_z_1)
w_x,w_z = grad(w,lu_x_1,piv_x_1,B_x_1,lu_z_1,piv_z_1,B_z_1)
b_x,b_z = grad(b,lu_x_1,piv_x_1,B_x_1,lu_z_1,piv_z_1,B_z_1)

g2u = grad2(u,lu_x_2,piv_x_2,B_x_2,lu_z_2,piv_z_2,B_z_2)
g2v = grad2(v,lu_x_2,piv_x_2,B_x_2,lu_z_2,piv_z_2,B_z_2)
g2w = grad2(w,lu_x_2,piv_x_2,B_x_2,lu_z_2,piv_z_2,B_z_2)
g2b = grad2(b,lu_x_2,piv_x_2,B_x_2,lu_z_2,piv_z_2,B_z_2)

rhs1 =  f*v  - (u*u_x + w*u_z) + nu*g2u
rhs2 = -f*u  - (u*v_x + w*v_z) + nu*g2v
rhs3 =    b  - (u*w_x + w*w_z) + nu*g2w
rhs4 = -N2*w - (u*b_x + w*b_z) + kappa*g2b



#---------------------------------------------------------------------------
# integrate the rhs terms one time step to get u*,b(tnp1) 
#---------------------------------------------------------------------------
ustar = u + dt*( rhs1 )
vstar = v + dt*( rhs2 )    
wstar = w + dt*( rhs3 )  
b     = b + dt*( rhs4 )

		   



#---------------------------------------------------------------------------
# determine the corrections needed for normal components at the boundaries 
# simply insert the BVALS for b 
# simply insert the BVALS for v  (only because d/dy=0 and v not pressure projected)
#---------------------------------------------------------------------------
#    EAST & WEST
for k in range(nz):
	U = parent_soln(x0,z0+z[k],tnp1,1)
	delta_u_e[k] = U - ustar[0,k]
	
	U = parent_soln(x0+Lx,z0+z[k],tnp1,1)
	delta_u_w[k] = U - ustar[nx-1,k]
	
	b[0,k]    = parent_soln(x0,z0+z[k],tnp1,4)
	b[nx-1,k] = parent_soln(x0+Lx,z0+z[k],tnp1,4)
	
	v[0,k]    = parent_soln(x0,z0+z[k],tnp1,2)
	v[nx-1,k] = parent_soln(x0+Lx,z0+z[k],tnp1,2)

#    BOTTOM & TOP	
for i in range(nx):
	W = parent_soln(x0+x[i],z0,tnp1,3)
	delta_w_b[i] = W - wstar[i,0]
	
	W = parent_soln(x0+x[i],z0+Lz,tnp1,3)
	delta_w_t[i] = W - wstar[i,nz-1]
	
	b[i,0]    = parent_soln(x0+x[i],z0,tnp1,4)
	b[i,nz-1] = parent_soln(x0+x[i],z0+Lz,tnp1,4)
	
	v[i,0]    = parent_soln(x0+x[i],z0,tnp1,2)
	v[i,nz-1] = parent_soln(x0+x[i],z0+Lz,tnp1,2)


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
	
	i=0
	ddx_v = ( -11.*vstar[i,k] + 18.*vstar[i+1,k] - 9.*vstar[i+2,k] + 2.* vstar[i+3,k] )/(6.*dx)
	ddn_v = -ddx_v
	delta_v_e[k] = gamma_x*( -V_x - ddn_v )
	
	i=0
	ddx_w = ( -11.*wstar[i,k] + 18.*wstar[i+1,k] - 9.*wstar[i+2,k] + 2.* wstar[i+3,k] )/(6.*dx)
	ddn_w = -ddx_w
	delta_w_e[k] = gamma_x*( -W_x - ddn_w )
	
	# WEST BOUNDARY  n_hat in positive x direction
	V_x,W_x = parent_derivs(x0+Lx,z0+z[k],tnp1,'x')
	
	i=nx-1
	ddx_v = ( 11.*vstar[i,k] - 18.*vstar[i-1,k] + 9.*vstar[i-2,k] - 2.* vstar[i-3,k] )/(6.*dx)
	ddn_v = ddx_v
	delta_v_w[k] = gamma_x*( V_x - ddn_v )
	
	i=nx-1
	ddx_w = ( 11.*wstar[i,k] - 18.*wstar[i-1,k] + 9.*wstar[i-2,k] - 2.* wstar[i-3,k] )/(6.*dx)
	ddn_w = ddx_w
	delta_w_w[k] = gamma_x*( W_x - ddn_w )

for i in range(nx):
	
	# BOTTOM BOUNDARY  n_hat in negative z direction
	U_z,V_z = parent_derivs(x0+x[i],z0,tnp1,'z')
	
	k=0
	ddz_u = ( -11.*ustar[i,k] + 18.*ustar[i,k+1] - 9.*ustar[i,k+2] + 2.* ustar[i,k+3] )/(6.*dz)
	ddn_u = -ddz_u
	delta_u_b[i] = gamma_z*( -U_z - ddn_u )
	
	ddz_v = ( -11.*vstar[i,k] + 18.*vstar[i,k+1] - 9.*vstar[i,k+2] + 2.* vstar[i,k+3] )/(6.*dz)
	ddn_v = -ddz_v
	delta_v_b[i] = gamma_z*( -V_z - ddn_v )
	
	# TOP BOUNDARY  n_hat in positive z direction
	U_z,V_z = parent_derivs(x0+x[i],z0+Lz,tnp1,'z')
	
	k=nz-1
	ddz_u = ( 11.*ustar[i,k] - 18.*ustar[i,k-1] + 9.*ustar[i,k-2] - 2.* ustar[i,k-3] )/(6.*dz)
	ddn_u = ddz_u
	delta_u_b[i] = gamma_z*( U_z - ddn_u )
	
	k=nz-1
	ddz_v = ( 11.*vstar[i,k] - 18.*vstar[i,k-1] + 9.*vstar[i,k-2] - 2.* vstar[i,k-3] )/(6.*dz)
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
div = divergence(ustar,wstar,lu_x_1,piv_x_1,B_x_1,lu_z_1,piv_z_1,B_z_1)
divmean = np.mean(div)
div = div - divmean


#---------------------------------------------------------------------------
#  invert the div of u_s^* to generate an approximate solution for phi
#  poisson equation has homogeneous Neumann conditions
#---------------------------------------------------------------------------
phi = poisson_invert_2d(div/dt,inversion_matrix)


#---------------------------------------------------------------------------
#  compute the gradient of phi
#---------------------------------------------------------------------------
phi_x,phi_z = grad(phi,lu_x_1,piv_x_1,B_x_1,lu_z_1,piv_z_1,B_z_1)


#---------------------------------------------------------------------------
#  do the pressure projection to get velocity components at tnp1
#---------------------------------------------------------------------------
u = ustar - dt*phi_x
w = wstar - dt*phi_z

		


#---------------------------------------------------------------------------
# compute exact soln at tnp1
#---------------------------------------------------------------------------
u_exact=np.zeros_like(u); v_exact=np.zeros_like(v) ; w_exact=np.zeros_like(w); b_exact=np.zeros_like(b)
for k in range(nz):
	for i in range(nx):
		
		#---------------------------------------------------------------------------
		# map nested to global coordinates, evaluate exact soln at tn  for use as IC 
		#---------------------------------------------------------------------------
		X = x0 + x[i]
		Z = z0 + z[k]	
		u_exact[i,k] = parent_soln(X,Z,tnp1,1)
		v_exact[i,k] = parent_soln(X,Z,tnp1,2)
		w_exact[i,k] = parent_soln(X,Z,tnp1,3)
		b_exact[i,k] = parent_soln(X,Z,tnp1,4)
				

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
fname = "pressure_test.nc"                                                                          # construct the netcdf filename
vars = [ustar,wstar,div,phi,u,v,w,b,u_exact,v_exact,w_exact,b_exact]                                # list of arrays with data
varnames = ['ustar','wstar','div','phi','u','v','w','b','u_exact','v_exact','w_exact','b_exact']    # list of variable names
dims = ['m/s','m/s','1/s','m2/s2','m/s','m/s','m/s','m/s2','m/s','m/s','m/s','m/s2']                # list of dimensions
success = write_netcdf_2d(vars,x,z,tnp1,fname,varnames,dims)              # create and write the netcdf file	
