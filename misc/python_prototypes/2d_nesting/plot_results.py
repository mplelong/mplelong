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

	
#w_profile = read_netcdf_2d_x0(varname,ncfile,tslice,nx-1)
#w_slice = read_netcdf_2d_z0(varname,ncfile,tslice,0)

tslice=0
slices = np.array([1024])

for idx in slices:
	ncfile = '1_period_output/XZ_' + str(idx).zfill(4) + '.nc'

	if( idx==slices[0] ):
		t,x,z = read_netcdf_2d_coords(ncfile)
		nx = x.size ; nz = z.size ; nt = t.size
		Lx = x[-1]-x[0]  ; Lz = z[-1]-z[0]
		I = (nx-1)/4   ;   K = (nz-1)/4
		
		fig = plt.figure(figsize=(8.5,11),dpi=600)    # fig size in inches
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
	
	wstar = read_netcdf_2d_x0('wstar',ncfile,tslice,I) ; wstar = wstar/np.max(np.abs(wstar))
	phi_z = read_netcdf_2d_x0('phi_z',ncfile,tslice,I) ; phi_z = phi_z/np.max(np.abs(phi_z))
	w = read_netcdf_2d_x0('w',ncfile,tslice,I)         ; w = w/np.max(np.abs(w))
	w_exact = read_netcdf_2d_x0('w_exact',ncfile,tslice,I)         ; w_exact = w_exact/np.max(np.abs(w_exact))
	
	ustar = read_netcdf_2d_x0('ustar',ncfile,tslice,I) ; ustar = ustar/np.max(np.abs(ustar))
	phi_x = read_netcdf_2d_x0('phi_x',ncfile,tslice,I) ; phi_x = phi_x/np.max(np.abs(phi_x))
	u = read_netcdf_2d_x0('u',ncfile,tslice,I)         ; u = u/np.max(np.abs(u))
	
	inc=8
	ax1 = plt.subplot(2,3,1)
	ax1.plot( wstar,z/Lz,'k',linewidth=2 ) #; ax1.plot( ustar,z/Lz,'k--' )
	ax2 = plt.subplot(2,3,2)
	ax2.plot( phi_z,z/Lz,'k',linewidth=2 ) #; ax2.plot( phi_x,z/Lz,'k--' )
	ax3 = plt.subplot(2,3,3)
	ax3.plot( w,z/Lz,'k',linewidth=2 )     #; ax3.plot( u,z/Lz,'k--' )
	ax3.plot( w_exact[0:nz:inc],z[0:nz:inc]/Lz,'ro',markersize=3 )
	
	
	wstar = read_netcdf_2d_z0('wstar',ncfile,tslice,K) ; wstar = wstar/np.max(np.abs(wstar))
	phi_z = read_netcdf_2d_z0('phi_z',ncfile,tslice,K) ; phi_z = phi_z/np.max(np.abs(phi_z))
	w = read_netcdf_2d_z0('w',ncfile,tslice,K)         ; w = w/np.max(np.abs(w))
	w_exact = read_netcdf_2d_z0('w_exact',ncfile,tslice,K)         ; w_exact = w_exact/np.max(np.abs(w_exact))
	
	ax1 = plt.subplot(6,1,4)
	ax1.plot( x/Lx,wstar,'k',linewidth=2 ) 
	ax2 = plt.subplot(6,1,5)
	ax2.plot( x/Lx,phi_z,'k',linewidth=2 )
	ax3 = plt.subplot(6,1,6)
	ax3.plot( x/Lx,w,'k',linewidth=2 )
	ax3.plot( x[0:nx:inc]/Lx,w_exact[0:nx:inc],'ro',markersize=3 )
	
	


#---------------------------------------------------------------------------------------
#   save the figure 
#--------------------------------------------------------------------------------------- 
fname = 'projection.eps' 
plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution

