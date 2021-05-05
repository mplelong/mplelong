#!
"""
    simple plot of cfl data, takes root_dir and tscale as a command line arguments
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from data_processing_utilities import parse_problem_params
#-----------------------------------------------------------------------------------------

root_dir = sys.argv[1]
root_dir = root_dir + '/'
tscale = sys.argv[2]

if(tscale=='s' or tscale=='secs'):
	xnorm = 1.
elif(tscale=='hrs'):
	xnorm = 3600.
elif(tscale=='days'):
	xnorm = (24.*3600.)
	
#--------------------------------------------------------------------------- 
# root directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
#root_dir = '/Users/kraig/flow_solve_BC/'
problem_params = parse_problem_params(root_dir)
[runlabel,restart_flag,do_second_scalar,AB_order,p1,p2,nx,ny,nz,dt,t0,tf,Lx,Ly,Lz, \
x_periodic,y_periodic,z_periodic,s1_def,s2_def,user_forcing_flag,rho0,g,f0,nu,     \
kappa1,kappa2,high_order_flag,p,T_diff] = problem_params

data_file = root_dir + 'output/cfl.dat'     # t, cfl_x, cfl_y, cfl_z
plot_dir  = root_dir + 'output/figures/'
plot_file = plot_dir + 'cfl.png'

if not os.path.exists(plot_dir):
	cmd = 'mkdir -p ' + plot_dir
	os.system(cmd)


#---------------------------------------
# Open file & read data columns...
#---------------------------------------
f = open(data_file, 'r')
f = np.loadtxt(data_file)
time = f[:,0] ; cfl_x = f[:,1]  ; cfl_y = f[:,2] ; cfl_z = f[:,3]




fig = plt.figure(figsize=(8,3.75),dpi=300)                 # fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])              # lower, bottom, width, height in ( 0 to 1)
axes.tick_params(direction='out', top=False, right=False,) # Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                      # Get rid of top axis line
axes.spines['right'].set_visible(False)

time = time/xnorm   # work in "tscale"

plt.plot(time,cfl_x,'r',label='cfl x')
plt.plot(time,cfl_y,'g',label='cfl y')
plt.plot(time,cfl_z,'b',label='cfl z')
plt.plot([time[0],time[-1]],[0.2,0.2],'k--')
plt.legend(loc='lower left')

axis_lims = [time[0], 1.05*time[-1], 0., 0.25]
plt.axis(axis_lims)

axes.set_xlabel(tscale,fontsize=14)
axes.set_ylabel(r'$(u,v,w)\,dt/(dx,dy,dz)$',fontsize=14)
title_string = "maximum cfl values vs time:          dt=%.4f  [s] " %(dt)
axes.set_title(title_string,fontsize=12)


plt.savefig(plot_file,dpi=300,bb_inches='tight')      # save plot file
