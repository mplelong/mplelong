#!python
"""
    Script to write global boundary data. Files created are
    south_vals.nc    south_derivs.nc
    north_vals.nc    north_derivs.nc
    east_vals.nc     east_derivs.nc
    west_vals.nc     west_derivs.nc
    bottom_vals.nc   bottom_derivs.nc
    top_vals.nc      top_derivs.nc
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,sys
import numpy as np
import scipy as sp
from scipy.io import netcdf  
from XZ_IWs_utilities import initialize_2d_variables,fill_boundary_vals,write_east_vals, \
                             write_west_vals,write_south_vals,write_north_vals,          \
                             write_bottom_vals,write_top_vals
#---------------------------------------------------------------------------------------

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
nu = 1.e-6                # [m2/s]        viscosity
kappa = nu                # [m2/s]        diffusivity
k_iw = 2.*pi/L            # [1/m]         horizontal wavenumber of internal wave
m_iw = 1.*pi/H            # [1/m]         vertical wavenumber of internal wave
omega2 = (k_iw**2 * N2 + m_iw**2 * f2)/(k_iw**2 + m_iw**2)
omega = np.sqrt(omega2)   # [1/s]         internal wave frequency
IWPER = (2.*pi/omega)     # [s]
A = 0.01                  # [m/s]         amplitude of U
phase = pi/4.             # [1]           initial phase shift of wave mode
WAVE_PARAMS = [A,f,N2,omega,k_iw,m_iw,phase]



#--------------------------------------------------------------------------- 
# get command line arguments
#---------------------------------------------------------------------------
top_dir = sys.argv[1]
#  add / for convenience
top_dir = top_dir + '/'

output_dir = top_dir + 'BVALS/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

nx = int(sys.argv[2])
ny = int(sys.argv[3])
nz = int(sys.argv[4])

nt = int(sys.argv[5])
dt = IWPER/(nt-1.)


print('...  create_global_bdry_files.py executing' )
print('...  flow_solve root directory:         ',top_dir )
print('...  global resolution nx x ny x nz:    ',nx,ny,nz )
print('...  number of time slices to save:     ',nt )
print('...  output directory for global boundary planes:   ',output_dir )

[U_bot,V_bot,W_bot,B_bot,U_top,V_top,W_top,B_top] = initialize_2d_variables(nx,ny)
[U_east,V_east,W_east,B_east,U_west,V_west,W_west,B_west] = initialize_2d_variables(ny,nz)
[U_south,V_south,W_south,B_south,U_north,V_north,W_north,B_north] = initialize_2d_variables(nx,nz)

[U_z_bot,V_z_bot,W_z_bot,B_z_bot,U_z_top,V_z_top,W_z_top,B_z_top] = initialize_2d_variables(nx,ny)
[U_x_east,V_x_east,W_x_east,B_x_east,U_x_west,V_x_west,W_x_west,B_x_west] = initialize_2d_variables(ny,nz)
[U_y_south,V_y_south,W_y_south,B_y_south,U_y_north,V_y_north,W_y_north,B_y_north] = initialize_2d_variables(nx,nz)

# pack things up for easy passing to functions
east_vals = [U_east,V_east,W_east,B_east]
west_vals = [U_west,V_west,W_west,B_west]
south_vals = [U_south,V_south,W_south,B_south]
north_vals = [U_north,V_north,W_north,B_north]
bot_vals =  [U_bot,V_bot,W_bot,B_bot]
top_vals =  [U_top,V_top,W_top,B_top]
BVALS = [east_vals,west_vals,south_vals,north_vals,bot_vals,top_vals]

east_derivs = [U_x_east,V_x_east,W_x_east,B_x_east] 
west_derivs = [U_x_west,V_x_west,W_x_west,B_x_west]
south_derivs = [U_y_south,V_y_south,W_y_south,B_y_south]
north_derivs = [U_y_north,V_y_north,W_y_north,B_y_north]
bot_derivs = [U_z_bot,V_z_bot,W_z_bot,B_z_bot]
top_derivs = [U_z_top,V_z_top,W_z_top,B_z_top] 
BDERIVS = [east_derivs,west_derivs,south_derivs,north_derivs,bot_derivs,top_derivs]

# north and south vals are not needed in this y periodic case, they can keep zero values

#---------------------------------------------------------------------------------------
#  parameters to define and control the nested simulation
#---------------------------------------------------------------------------------------
Lx = 30000.                      # [m]      size of nested domain in x 
Lz = 600.                        # [m]      size of nested domain in z
Ly = Lx
x0 = 0.50*L				         # [m]      horizontal offset for origin of nested domain
y0 = 0.40*L				         # [m]      horizontal offset for origin of nested domain
z0 = 0.60*H                      # [m]      vertical offset for origin of nested domain

# x nested grid 	
dx = Lx/(nx-1.)                  # [m] grid spacing
x = np.linspace(0,Lx,nx)         # [m] discrete grid points in [0,Lx]

# y nested grid    	
dy = Ly/(ny-1.)                  # [m] grid spacing
y = np.linspace(0,Ly,ny)         # [m] discrete grid points in [0,Lx]

# z nested grid 	
dz = Lz/(nz-1.)                  # [m] grid spacing
z = np.linspace(0,Lz,nz)         # [m] discrete grid points in [0,Lz]


for islice in range(nt):

	t = islice*dt   # current time in [s]
	
	#------------------------------------------------------------------------------------
	#  Get the boundary values at time t and store in BVALS and BDERIVS. 
	#------------------------------------------------------------------------------------
	BVALS,BDERIVS = fill_boundary_vals(x,y,z,t,x0,y0,z0,WAVE_PARAMS,BVALS,BDERIVS)

	# unpack boundary information into 1d arrays
	[east_vals,west_vals,south_vals,north_vals,bot_vals,top_vals]=BVALS
	[U_east,V_east,W_east,B_east] = east_vals
	[U_west,V_west,W_west,B_west] = west_vals
	[U_south,V_south,W_south,B_south] = south_vals
	[U_north,V_north,W_north,B_north] = north_vals
	[U_bot,V_bot,W_bot,B_bot] = bot_vals
	[U_top,V_top,W_top,B_top] = top_vals
	
	[east_derivs,west_derivs,south_derivs,north_derivs,bot_derivs,top_derivs] = BDERIVS
	[U_x_east,V_x_east,W_x_east,B_x_east] = east_derivs
	[U_x_west,V_x_west,W_x_west,B_x_west] = west_derivs
	[U_y_south,V_y_south,W_y_south,B_y_south] = south_derivs
	[U_y_north,V_y_north,W_y_north,B_y_north] = north_derivs
	[U_z_bot,V_z_bot,W_z_bot,B_z_bot] = bot_derivs
	[U_z_top,V_z_top,W_z_top,B_z_top] = top_derivs
	
	success = write_east_vals(t,islice,y,z,east_vals,east_derivs,output_dir,nt)
	success = write_west_vals(t,islice,y,z,west_vals,west_derivs,output_dir,nt)
	
	success = write_south_vals(t,islice,x,z,south_vals,south_derivs,output_dir,nt)
	success = write_north_vals(t,islice,x,z,north_vals,north_derivs,output_dir,nt)
	
	success = write_bottom_vals(t,islice,x,y,bot_vals,bot_derivs,output_dir,nt)
	success = write_top_vals(t,islice,x,y,top_vals,top_derivs,output_dir,nt)
	
	if(np.mod(islice,5)==0):
		print('...   wrote E/W, S/N and B/T boundary data at time slice ',islice )

command = 'ls -lh ' + output_dir
os.system(command) 
exit()
