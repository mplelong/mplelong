"""
Script to read in 2 consecutive time slices from a global 2d netcdf file and make
a sequence of images showing the progression through the projection algorithm
         w -> w*^ -> div -> phi -> phi_z -> w
"""
import os,math,sys              
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from netcdf_stuff import read_netcdf_2d_coords,read_netcdf_2d_var
from data_processing_utilities import prepare_BernoulliCosine,evaluate_basis_functions,grad_2d

data_dir = '../../output/2D/'
plot_dir = './'
do_print_eps = True
do_print_png = True


#---------------------------------------------------------------------------------------
#  parameters for the environment and the outer solution which is used to
#  supply initial and time dependent boundary conditions
#---------------------------------------------------------------------------------------
pi = np.pi
H = 3000.                 # [m]           full water depth
L = 1.5e5                 # [m]           horizontal scale = 150 km
N = 2.0e-3                # [1/s]         buoyancy frequency 
BVPER = (2.*pi/N)/3600.   # [hours]       buoyancy period
f = 1.e-4                 # [1/s]         Coriolis parameter
IPER = (2.*pi/f)/3600.    # [hours]       inertial period
N2 = N*N; f2 = f*f        # [1/s2]        squared frequencies

k_iw = 2.*pi/L            # [1/m]         horizontal wavenumber of internal wave
m_iw = 1.*pi/H            # [1/m]         vertical wavenumber of internal wave

omega2 = (k_iw**2 * N2 + m_iw**2 * f2)/(k_iw**2 + m_iw**2)
omega = np.sqrt(omega2)   # [1/s]         internal wave frequency
IWPER = (2.*pi/omega)     # [s]
A = 0.01                  # [m/s]         amplitude of U
phase = pi/4.             # [1]           initial phase shift of wave mode


def parent_soln(X,Z,t,id):
	#-----------------------------------------------------------------
	# function to evaluate parent soln at X,Z,t  ; 
	# here an IW mode with wave parameters available globally
	# id = [0,1,2,3,6] for [U,V,W,B,P] respectively
	#-----------------------------------------------------------------
	argz = m_iw*Z
	argx = k_iw*X - omega*t + phase
	if(id==0): ans = A*np.cos(argz)*np.cos(argx) # U
	if(id==1): ans = A*(f/omega)*np.cos(argz)*np.sin(argx) # V    
	if(id==2): ans = A*(k_iw/m_iw)*np.sin(argz)*np.sin(argx) # W 
	if(id==3): ans = -A*(k_iw/m_iw)*(N2/omega)*np.sin(argz)*np.cos(argx) # B
	if(id==6): ans = -A*(1./(k_iw*omega))*(f2-omega2)*np.cos(argz)*np.cos(argx) # P	
	return ans





filename = data_dir + 'XZ_0.nc'
#------------------------------------------------------------------
#  get the nested domain coordinate data
#------------------------------------------------------------------
x0 = 0.50*L				                         # [m]  horizontal offset for origin of nested domain
z0 = 0.60*H                                      # [m]  vertical offset for origin of nested domain
t,x,z = read_netcdf_2d_coords(filename,'x','z')  # [s],[m],[m]
[nx,nz,nt] = [x.size,z.size,t.size]              # nx, nz, nt
Lx = x[-1]-x[0] ; dx = x[1]-x[0]                 # Lx  dx
Lz = z[-1]-z[0] ; dz = z[1]-z[0]                 # Lz  dz    
DT = t[1]-t[0]                                   # time increment between saved slices

#--------------------------------------------------------------------------
# evaluate the basis functions for the S_0 and S_L (U_n/B_n+1) series
# and for their derivatives 
#   x_basis_functions=[S_0,S_L,S_0_x,S_L_x] similar for z
#--------------------------------------------------------------------------
Q = 9
frac = 0.075
x_basis_functions = evaluate_basis_functions(x,Q)
z_basis_functions = evaluate_basis_functions(z,Q)
basis_functions = [x_basis_functions,z_basis_functions]

#--------------------------------------------------------------------------
# create and factor the matrices for solving for the expansion coeffs
# in the Bernoulli/Cosine differentiation scheme
#--------------------------------------------------------------------------
LU_x = prepare_BernoulliCosine(x,Q)
LU_z = prepare_BernoulliCosine(z,Q)
print('............................................     BernoulliCosine prep done ')




for islice in np.array([nt-2]):
	varname = 'w'
	w = read_netcdf_2d_var(varname,filename,islice)        # w[z,x]
	
	varname = 'wstar'
	wstar = read_netcdf_2d_var(varname,filename,islice)    # wstar[z,x]
	
	varname = 'divustar'
	divustar = read_netcdf_2d_var(varname,filename,islice) # div_u*[z,x]
	
	varname = 'phi'
	phi = read_netcdf_2d_var(varname,filename,islice)      # phi[z,x]
	
	varname = 'w'
	wnp1 = read_netcdf_2d_var(varname,filename,islice+1)   # w[z,x]  at tn + dt
	
	# compute pressure gradient;   grad_2d takes 2d input f[x,z]
	[phi_x,phi_z] = grad_2d(phi.transpose(),x,z,frac,LU_x,LU_z,Q,basis_functions)
	phi_x = phi_x.transpose()
	phi_z = phi_z.transpose()
	print('............................................     pressure gradient computed ')
	

	#---------------------------------------------------------------------------------------
	#  generic setup for the figure
	#---------------------------------------------------------------------------------------
	FS = 10                                                   # font size
	fig = plt.figure(figsize=(8,4),dpi=600)                 # fig size in inches
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	axes = fig.add_axes([0.10, 0.15, 0.85, 0.85])               # 
	axes.set_facecolor('xkcd:light grey')
	XmajorLocator   = MultipleLocator(10)            # major ticks x
	XmajorFormatter = FormatStrFormatter('%d')       # %d is integer  %.2f
	XminorLocator   = MultipleLocator(5)             # minor ticks x

	YmajorLocator   = MultipleLocator(0.1)
	YmajorFormatter = FormatStrFormatter('%.2f')
	YminorLocator   = MultipleLocator(0.25)
	axis_lims = [x0/1000., (x0+Lx)/1000., z0/1000.,(z0+Lz)/1000.]
	
	
	#---------------------------------------------------------------------------------------
	# w(tn)
	#---------------------------------------------------------------------------------------
	
	ax1 = plt.subplot(2,3,1)
	#---------------------------------------------------------------------------------------
	# axis limits for plot and tick marks
	#---------------------------------------------------------------------------------------
	ax1.tick_params(direction='out', top=False, right=False,)  # Turn ticks out
	ax1.tick_params(axis='both', which='major', labelsize=FS)
	ax1.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
	ax1.spines['right'].set_visible(True)
	
	v_min=np.min(w) ; v_max=-v_min
	CMAP = plt.get_cmap('bwr')    # bwr     
	ncontours = 128 ; levels = np.linspace(v_min,v_max,ncontours)
	
	#  plot w at start of time step
	plt.contour((x0+x)/1000.,(z0+z)/1000.,w,levels,cmap=CMAP,extend='both',zorder=2,linewidths=0.5) 
	plt.clim(v_min,v_max)

               
	#---------------------------------------------------------------------------------------
	# add tick marks, impose axis limits
	#---------------------------------------------------------------------------------------
	XmajorFormatter = FormatStrFormatter('')    # %d
	ax1.axis(axis_lims)
	ax1.xaxis.set_major_locator(XmajorLocator)
	ax1.xaxis.set_major_formatter(XmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.xaxis.set_minor_locator(XminorLocator)

	ax1.yaxis.set_major_locator(YmajorLocator)
	ax1.yaxis.set_major_formatter(YmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.yaxis.set_minor_locator(YminorLocator)

	#---------------------------------------------------------------------------------------
	#  axis labels and title
	#---------------------------------------------------------------------------------------
	#ax1.set_xlabel(r'$X \, \rm{[km]}$',fontsize=FS+2)
	ax1.set_ylabel(r'$Z \, \rm{[km]}$',fontsize=FS+2)
	title_string = r"$w(t^n)$" 
	ax1.set_title(title_string, fontsize=FS+2)
	
	
	#---------------------------------------------------------------------------------------
	# w*hat(tn)
	#---------------------------------------------------------------------------------------
	
	ax1 = plt.subplot(2,3,2)
	#---------------------------------------------------------------------------------------
	# axis limits for plot and tick marks
	#---------------------------------------------------------------------------------------
	ax1.tick_params(direction='out', top=False, right=False,)  # Turn ticks out
	ax1.tick_params(axis='both', which='major', labelsize=FS)
	ax1.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
	ax1.spines['right'].set_visible(True)
	
	v_max=np.max(np.abs(wstar)) ; v_min=-v_max 
	CMAP = plt.get_cmap('bwr')    # bwr     
	ncontours = 64 ; levels = np.linspace(v_min,v_max,ncontours)
	
	#  plot w*hat
	plt.contour((x0+x)/1000.,(z0+z)/1000.,wstar,levels,cmap=CMAP,extend='both',zorder=2,linewidths=0.5) 
	plt.clim(v_min,v_max)

               
	#---------------------------------------------------------------------------------------
	# add tick marks, impose axis limits
	#---------------------------------------------------------------------------------------
	XmajorFormatter = FormatStrFormatter('')    # %d
	ax1.axis(axis_lims)
	ax1.xaxis.set_major_locator(XmajorLocator)
	ax1.xaxis.set_major_formatter(XmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.xaxis.set_minor_locator(XminorLocator)

	YmajorFormatter = FormatStrFormatter('')
	ax1.yaxis.set_major_locator(YmajorLocator)
	ax1.yaxis.set_major_formatter(YmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.yaxis.set_minor_locator(YminorLocator)

	#---------------------------------------------------------------------------------------
	#  axis labels and title
	#---------------------------------------------------------------------------------------
	#ax1.set_xlabel(r'$X \, \rm{[km]}$',fontsize=FS+2)
	#ax1.set_ylabel(r'$Z \, \rm{[m]}$',fontsize=FS+2)
	title_string = r"${\hat w}_*$" 
	ax1.set_title(title_string, fontsize=FS+2)
	
	
	
	
	#---------------------------------------------------------------------------------------
	# divustar
	#---------------------------------------------------------------------------------------
	
	ax1 = plt.subplot(2,3,3)
	#---------------------------------------------------------------------------------------
	# axis limits for plot and tick marks
	#---------------------------------------------------------------------------------------
	ax1.tick_params(direction='out', top=False, right=False,)  # Turn ticks out
	ax1.tick_params(axis='both', which='major', labelsize=FS)
	ax1.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
	ax1.spines['right'].set_visible(True)
	
	v_max=np.max(np.abs(divustar)) ; v_min=-v_max 
	CMAP = plt.get_cmap('bwr')    # bwr     
	ncontours = 16 ; levels = np.linspace(v_min,v_max,ncontours)
	
	#  plot divustar
	plt.contour((x0+x)/1000.,(z0+z)/1000.,divustar,levels,cmap=CMAP,extend='both',zorder=2,linewidths=0.5) 
	plt.clim(v_min,v_max)

               
	#---------------------------------------------------------------------------------------
	# add tick marks, impose axis limits
	#---------------------------------------------------------------------------------------
	XmajorFormatter = FormatStrFormatter('')    # %d
	ax1.axis(axis_lims)
	ax1.xaxis.set_major_locator(XmajorLocator)
	ax1.xaxis.set_major_formatter(XmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.xaxis.set_minor_locator(XminorLocator)

	YmajorFormatter = FormatStrFormatter('')
	ax1.yaxis.set_major_locator(YmajorLocator)
	ax1.yaxis.set_major_formatter(YmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.yaxis.set_minor_locator(YminorLocator)

	#---------------------------------------------------------------------------------------
	#  axis labels and title
	#---------------------------------------------------------------------------------------
	#ax1.set_xlabel(r'$X \, \rm{[km]}$',fontsize=FS+2)
	#ax1.set_ylabel(r'$Z \, \rm{[m]}$',fontsize=FS+2)
	title_string = r"$\nabla \cdot {\vec {\hat u}}_*$" 
	ax1.set_title(title_string, fontsize=FS+2)
	
	
	
	#---------------------------------------------------------------------------------------
	# phi
	#---------------------------------------------------------------------------------------	
	ax1 = plt.subplot(2,3,4)
	#---------------------------------------------------------------------------------------
	# axis limits for plot and tick marks
	#---------------------------------------------------------------------------------------
	ax1.tick_params(direction='out', top=False, right=False,)  # Turn ticks out
	ax1.tick_params(axis='both', which='major', labelsize=FS)
	ax1.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
	ax1.spines['right'].set_visible(True)
	
	v_min=np.min(phi) ; v_max=np.max(phi)
	CMAP = plt.get_cmap('bwr')    # bwr     
	ncontours = 64 ; levels = np.linspace(v_min,v_max,ncontours)
	
	#  plot phi
	plt.contour((x0+x)/1000.,(z0+z)/1000.,phi,levels,cmap=CMAP,extend='both',zorder=2,linewidths=0.5) 
	plt.clim(v_min,v_max)

               
	#---------------------------------------------------------------------------------------
	# add tick marks, impose axis limits
	#---------------------------------------------------------------------------------------
	XmajorFormatter = FormatStrFormatter('%d')
	ax1.axis(axis_lims)
	ax1.xaxis.set_major_locator(XmajorLocator)
	ax1.xaxis.set_major_formatter(XmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.xaxis.set_minor_locator(XminorLocator)

	YmajorFormatter = FormatStrFormatter('%.2f')
	ax1.yaxis.set_major_locator(YmajorLocator)
	ax1.yaxis.set_major_formatter(YmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.yaxis.set_minor_locator(YminorLocator)

	#---------------------------------------------------------------------------------------
	#  axis labels and title
	#---------------------------------------------------------------------------------------
	ax1.set_xlabel(r'$X \, \rm{[km]}$',fontsize=FS+2)
	ax1.set_ylabel(r'$Z \, \rm{[km]}$',fontsize=FS+2)
	title_string = r"$\phi$" 
	ax1.set_title(title_string, fontsize=FS+2)
	
	
	
	#---------------------------------------------------------------------------------------
	# phi_z
	#---------------------------------------------------------------------------------------	
	ax1 = plt.subplot(2,3,5)
	#---------------------------------------------------------------------------------------
	# axis limits for plot and tick marks
	#---------------------------------------------------------------------------------------
	ax1.tick_params(direction='out', top=False, right=False,)  # Turn ticks out
	ax1.tick_params(axis='both', which='major', labelsize=FS)
	ax1.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
	ax1.spines['right'].set_visible(True)
	
	v_max=np.max(np.abs(phi_z)) ; v_min=-v_max 
	CMAP = plt.get_cmap('bwr')    # bwr     
	ncontours = 64 ; levels = np.linspace(v_min,v_max,ncontours)
	
	#  plot phi
	plt.contour((x0+x)/1000.,(z0+z)/1000.,phi_z,levels,cmap=CMAP,extend='both',zorder=2,linewidths=0.5) 
	plt.clim(v_min,v_max)

               
	#---------------------------------------------------------------------------------------
	# add tick marks, impose axis limits
	#---------------------------------------------------------------------------------------
	XmajorFormatter = FormatStrFormatter('%d')
	ax1.axis(axis_lims)
	ax1.xaxis.set_major_locator(XmajorLocator)
	ax1.xaxis.set_major_formatter(XmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.xaxis.set_minor_locator(XminorLocator)

	YmajorFormatter = FormatStrFormatter('')
	ax1.yaxis.set_major_locator(YmajorLocator)
	ax1.yaxis.set_major_formatter(YmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.yaxis.set_minor_locator(YminorLocator)

	#---------------------------------------------------------------------------------------
	#  axis labels and title
	#---------------------------------------------------------------------------------------
	ax1.set_xlabel(r'$X \, \rm{[km]}$',fontsize=FS+2)
	#ax1.set_ylabel(r'$Z \, \rm{[km]}$',fontsize=FS+2)
	title_string = r"$\phi_z$" 
	ax1.set_title(title_string, fontsize=FS+2)
	
	
	
	
	
	
	
	
	#---------------------------------------------------------------------------------------
	# w at tnp1
	#---------------------------------------------------------------------------------------
	
	ax1 = plt.subplot(2,3,6)
	#---------------------------------------------------------------------------------------
	# axis limits for plot and tick marks
	#---------------------------------------------------------------------------------------
	ax1.tick_params(direction='out', top=False, right=False,)  # Turn ticks out
	ax1.tick_params(axis='both', which='major', labelsize=FS)
	ax1.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
	ax1.spines['right'].set_visible(True)
	
	v_min=np.min(w) ; v_max=-v_min
	CMAP = plt.get_cmap('bwr')    # bwr     
	ncontours = 128 ; levels = np.linspace(v_min,v_max,ncontours)
	
	#  plot phi
	plt.contour((x0+x)/1000.,(z0+z)/1000.,wnp1,levels,cmap=CMAP,extend='both',zorder=2,linewidths=0.5) 
	plt.clim(v_min,v_max)

               
	#---------------------------------------------------------------------------------------
	# add tick marks, impose axis limits
	#---------------------------------------------------------------------------------------
	ax1.axis(axis_lims)
	ax1.xaxis.set_major_locator(XmajorLocator)
	ax1.xaxis.set_major_formatter(XmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.xaxis.set_minor_locator(XminorLocator)

	YmajorFormatter = FormatStrFormatter('')
	ax1.yaxis.set_major_locator(YmajorLocator)
	ax1.yaxis.set_major_formatter(YmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	ax1.yaxis.set_minor_locator(YminorLocator)

	#---------------------------------------------------------------------------------------
	#  axis labels and title
	#---------------------------------------------------------------------------------------
	ax1.set_xlabel(r'$X \, \rm{[km]}$',fontsize=FS+2)
	#ax1.set_ylabel(r'$Z \, \rm{[km]}$',fontsize=FS+2)
	title_string = r"$w(t^{n+1})$" 
	ax1.set_title(title_string, fontsize=FS+2)
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	#---------------------------------------------------------------------------------------
	#  possibly, print the figure as an eps file
	#---------------------------------------------------------------------------------------
	if do_print_eps == 1:
	    fname = plot_dir + 'w_projection.eps'
	    plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution

	#---------------------------------------------------------------------------------------
	#  possibly, print the figure as a png file
	#---------------------------------------------------------------------------------------
	if do_print_png == 1:
	    fname = plot_dir + 'w_projection.png'
	    plt.savefig(fname,dpi=600,bb_inches='tight') # save a png
	    
	print('... saved the image file ',fname )


	plt.close()
	
exit()


