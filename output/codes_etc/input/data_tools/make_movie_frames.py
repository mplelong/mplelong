"""
Script to read in time slices from a global 2d netcdf file and make
plots of the flow_solve solution embedded within the outer domain
"""
import os,math,sys              
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from netcdf_stuff import read_netcdf_2d_coords,read_netcdf_2d_var

data_dir = '../../output/2D/'
plot_dir = '../../output/movie_frames/'
do_print_eps = False
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

NX=64 ; DX= L/NX; X=np.linspace(0,L-DX,NX)     # open interval   [0,L)
NZ=65 ; Z = np.linspace(0,H,NZ)                # closed interval [0,H]


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

for islice in range(nt):
	frame = islice
	varname = 'u'
	u = read_netcdf_2d_var(varname,filename,islice)    # u[z,x]

	U = np.zeros( (NZ,NX), float )
	for i in range(NX):
		for k in range(NZ):
			varid = 0  # U
			U[k,i] = parent_soln(X[i],Z[k],t[islice],varid)


	#---------------------------------------------------------------------------------------
	#  generic setup for the figure
	#---------------------------------------------------------------------------------------
	FS = 12                                                   # font size
	fig = plt.figure(figsize=(8,4),dpi=600)                 # fig size in inches
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	axes = fig.add_axes([0.10, 0.20, 0.85, 0.70])               # 
	axes.tick_params(direction='out', top=False, right=False,)  # Turn ticks out
	axes.tick_params(axis='both', which='major', labelsize=FS)
	axes.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
	axes.spines['right'].set_visible(False)
	axes.set_facecolor('xkcd:light grey')

	#---------------------------------------------------------------------------------------
	# define domain extent for image map, axis limits for plot and tick marks
	#---------------------------------------------------------------------------------------
	axis_lims = [0., L/1000., H/2, H]

	XmajorLocator   = MultipleLocator(25)            # major ticks x
	XmajorFormatter = FormatStrFormatter('%d')       # %d is integer  %.2f
	XminorLocator   = MultipleLocator(5)             # minor ticks x

	YmajorLocator   = MultipleLocator(500)
	YmajorFormatter = FormatStrFormatter('%d')
	YminorLocator   = MultipleLocator(250)

	v_min=-0.01; v_max=0.01
	CMAP = plt.get_cmap('bwr')    # bwr     
	ncontours = 32 ; levels = np.linspace(v_min,v_max,ncontours)
	plt.contour(X/1000.,Z,U,levels,cmap=CMAP,extend='both',zorder=1,linewidths=0.75) 
	plt.clim(v_min,v_max)
    
	#---------------------------------------------------------------------------------------
	#  plot the soln in the nested domain
	#---------------------------------------------------------------------------------------
	plt.plot([x0/1000.,(x0+Lx)/1000.],[z0,z0],'k')
	plt.plot([x0/1000.,(x0+Lx)/1000.],[z0+Lz,z0+Lz],'k')
	plt.plot([x0/1000., x0/1000.],[z0,z0+Lz],'k')
	plt.plot([(x0+Lx)/1000., (x0+Lx)/1000.],[z0,z0+Lz],'k')
	plt.contour((x0+x)/1000.,(z0+z),u,levels,cmap=CMAP,extend='both',zorder=2,linewidths=3.0) 
	plt.clim(v_min,v_max)

	plt.contour((x0+x)/1000.,(z0+z),u,levels,colors='k',extend='both',zorder=3,linewidths=0.5) 
               
	#---------------------------------------------------------------------------------------
	# add tick marks, impose axis limits
	#---------------------------------------------------------------------------------------
	plt.axis(axis_lims)

	axes.xaxis.set_major_locator(XmajorLocator)
	axes.xaxis.set_major_formatter(XmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	axes.xaxis.set_minor_locator(XminorLocator)

	axes.yaxis.set_major_locator(YmajorLocator)
	axes.yaxis.set_major_formatter(YmajorFormatter)
	#for the minor ticks, use no labels; default NullFormatter
	axes.yaxis.set_minor_locator(YminorLocator)

	#---------------------------------------------------------------------------------------
	#  axis labels and title
	#---------------------------------------------------------------------------------------
	axes.set_xlabel(r'$X \, \rm{[km]}$',fontsize=FS+2)
	axes.set_ylabel(r'$Z \, \rm{[m]}$',fontsize=FS+2)
	title_string = r"$u(x,z,t)$   t=%.2f IW periods " %(t[islice]/IWPER)
	axes.set_title(title_string, fontsize=FS+2)

	#---------------------------------------------------------------------------------------
	#  possibly, print the figure as an eps file
	#---------------------------------------------------------------------------------------
	if do_print_eps == 1:
	    fname = plot_dir + 'U_' + str(frame).zfill(4) + '.eps'
	    plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution

	#---------------------------------------------------------------------------------------
	#  possibly, print the figure as a png file
	#---------------------------------------------------------------------------------------
	if do_print_png == 1:
	    fname = plot_dir + 'U_' + str(frame).zfill(4) + '.png'
	    plt.savefig(fname,dpi=600,bb_inches='tight') # save a png
	    
	print('... saved the image file ',fname )


	plt.close()
	
exit()


