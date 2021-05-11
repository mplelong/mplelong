#!
"""
    simple plot to compare w/ exact solution, takes root_dir as a command line argument
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
	
#--------------------------------------------------------------------------- 
# root directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
#root_dir = '/Users/kraig/flow_solve_BC/'
problem_params = parse_problem_params(root_dir)
[runlabel,restart_flag,do_second_scalar,AB_order,p1,p2,nx,ny,nz,dt,t0,tf,Lx,Ly,Lz, \
x_periodic,y_periodic,z_periodic,s1_def,s2_def,user_forcing_flag,rho0,g,f0,nu,     \
kappa1,kappa2,high_order_flag,p,T_diff] = problem_params

data_file = root_dir + 'output/TG_results'     # xval,yval,zval,t(n),uvw(n,1),uvw(n,2),uvw(n,3),exact(n,1),exact(n,2),exact(n,3)
plot_dir  = root_dir + 'output/figures/'
plot_file = plot_dir + 'comparison.pdf'

if not os.path.exists(plot_dir):
	cmd = 'mkdir -p ' + plot_dir
	os.system(cmd)


#---------------------------------------
# Open file & read data columns...
#---------------------------------------
f = open(data_file, 'r')
f = np.loadtxt(data_file)
time = f[:,3] ; u = f[:,4]  ; v = f[:,5] ; w = f[:,6]
u_exact = f[:,7]  ; v_exact = f[:,8] ; w_exact = f[:,9]
nt=time.size


# Update the matplotlib configuration parameters:
FS=10
matplotlib.rcParams.update({'font.size': FS, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})

fig = plt.figure(figsize=(6,3),dpi=150)                 # fig size in inches
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

axes = fig.add_axes([0.15, 0.15, 0.75, 0.75])              # lower, bottom, width, height in ( 0 to 1)
axes.tick_params(direction='out', top=False, right=False,) # Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                      # Get rid of top axis line
axes.spines['right'].set_visible(False)


inc=25; MS=2
plt.plot(time[0:nt:inc],u[0:nt:inc],'ro',time,u_exact,'r',markersize=MS)
plt.plot(time[0:nt:inc],v[0:nt:inc],'go',time,v_exact,'g',markersize=MS)
plt.plot(time[0:nt:inc],w[0:nt:inc],'bo',time,w_exact,'b',markersize=MS)


axis_lims = [time[0], 1.05*time[-1], 0., 1.05]
#plt.axis(axis_lims)

axes.set_xlabel('time [s]',fontsize=FS)
axes.set_ylabel(r'[m/s]',fontsize=FS)
title_string = "computed vs exact soln at arbitrary point: 128x128x17  L=1 [m] nu=1.d-2 [m2/s]  dt=%.4f  [s] " %(dt)
axes.set_title(title_string,fontsize=FS-2)


plt.savefig(plot_file,dpi=300)      # save plot file
