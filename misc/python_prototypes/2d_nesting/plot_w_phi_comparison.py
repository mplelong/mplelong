#!/usr/local/bin/python
"""
  Illustrate/explore the projection scheme for boundary forced problems
  
  Read and plot the netcdf data saved by projection_step_details.py
  
  1d arrays
  	t,x,z    (these files have only 1 time slice)
  2d arrays 
  	netcdf storage convention:  u[t,z,x] etc  (different than projection_step_details.py)
  	u,v,w,b,ustar,vstar,wstar,u_exact,v_exact,w_exact,b_exact,phi,phi_x,phi_z
      
"""
import os,math,sys              
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
	
from netcdf_stuff import read_netcdf_2d_var, read_netcdf_2d_coords
from netcdf_stuff import read_netcdf_2d_x0, read_netcdf_2d_z0

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

x0 = 0.50*L				         # [m]      horizontal offset for origin of nested domain
z0 = 0.60*H                      # [m]      vertical offset for origin of nested domain


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
	
	
	
	
	
	
	
	

#w_profile = read_netcdf_2d_x0(varname,ncfile,tslice,nx-1)
#w_slice = read_netcdf_2d_z0(varname,ncfile,tslice,0)

tslice=0
w_plot = 1  ; p_plot = 2
slices = np.array([256,512,768,1024])
for idx in slices:
	ncfile = '1_period_output/XZ_' + str(idx).zfill(4) + '.nc'

	if( idx==slices[0] ):
		t,x,z = read_netcdf_2d_coords(ncfile)
		nx = x.size ; nz = z.size ; nt = t.size
		Lx = x[-1]-x[0]  ; Lz = z[-1]-z[0]
		
		fig = plt.figure(figsize=(8.5,11),dpi=600)    # fig size in inches
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
	
	w = read_netcdf_2d_var('w',ncfile,tslice)
	w_exact = read_netcdf_2d_var('w_exact',ncfile,tslice)
	
	phi = read_netcdf_2d_var('phi',ncfile,tslice)
	
	phi_exact = np.zeros_like(phi)          # forgot to save so just recompute
	t,x,z = read_netcdf_2d_coords(ncfile)   # need t
	for k in range(nz):
		for i in range(nx):
			X = x0+x[i] ; Z=z0+z[k]
			phi_exact[k,i] = parent_soln(X,Z,t,6)   # note switch in indices
	
	phi = phi - np.mean(phi)
	phi_exact = phi_exact - np.mean(phi_exact)
	
	if( idx==slices[0] ): 
		W0 = np.max( np.abs(w_exact) )
		P0 = np.max( np.abs(phi) )



	axw = plt.subplot(4,2,w_plot)
	nlevels=256
	axw.contourf(x/Lx, z/Lz, w/W0, nlevels, cmap = 'bwr' )   
	nlevels=16
	axw.contour(x/Lx, z/Lz, w/W0, nlevels, colors='black', linewidths=2 )
	axw.contour(x/Lx, z/Lz, w_exact/W0, nlevels, colors='yellow', linestyles='dotted',linewidths=2)
	w_plot = w_plot+2
	
	axw = plt.subplot(4,2,p_plot)
	nlevels=256
	axw.contourf(x/Lx, z/Lz, phi/P0, nlevels, cmap = 'PiYG' )   
	nlevels=16
	axw.contour(x/Lx, z/Lz, phi/P0, nlevels, colors='black', linewidths=2 )
	axw.contour(x/Lx, z/Lz, phi_exact/P0, nlevels, colors='yellow', linestyles='dotted',linewidths=2)
	p_plot = p_plot+2


#---------------------------------------------------------------------------------------
#   save the figure 
#--------------------------------------------------------------------------------------- 
fname = 'w_phi_comparison.eps' 
plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution

