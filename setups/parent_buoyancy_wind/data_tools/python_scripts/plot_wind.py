#!
"""
    simple plot of wind speed and direction, takes root_dir 
    and tscale as a command line arguments
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
print(root_dir)
tscale = sys.argv[2]
print(tscale)

if(tscale=='s'):
	xnorm = 1.
elif(tscale=='hrs'):
	xnorm = 3600.
elif(tscale=='days'):
	xnorm = (24.*3600.)
	
#--------------------------------------------------------------------------- 
# root directory 
# (below which the output directory containing the saved data files exists)
#---------------------------------------------------------------------------
problem_params = parse_problem_params(root_dir)
[runlabel,restart_flag,do_second_scalar,AB_order,p1,p2,nx,ny,nz,dt,t0,tf,Lx,Ly,Lz, \
x_periodic,y_periodic,z_periodic,s1_def,s2_def,user_forcing_flag,rho0,g,f0,nu,     \
kappa1,kappa2,high_order_flag,p,T_diff] = problem_params

data_file = root_dir + 'output/1D/wind_speed_direction'     # t [s], u_10 [m/s], theta [rad]
plot_dir  = root_dir + 'output/figures/'
plot_file = plot_dir + 'wind.pdf'

if not os.path.exists(plot_dir):
	cmd = 'mkdir -p ' + plot_dir
	os.system(cmd)


#---------------------------------------
# Open file & read data columns...
#---------------------------------------
f = open(data_file, 'r')
f = np.loadtxt(data_file)
time = f[:,0] ; speed = f[:,1]  ; dir = f[:,2] 



FS=14
fig = plt.figure(figsize=(8,3.75),dpi=300)                 # fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])              # lower, bottom, width, height in ( 0 to 1)
axes.tick_params(direction='out', top=False, right=False,) # Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                      # Get rid of top axis line
axes.spines['right'].set_visible(False)

time = time/xnorm   # work in "tscale"
mean = 12.0         # prescribed mean wind speed


LC='r'
axes.set_xlabel(r'$t \, \rm[days]$',fontsize=FS)
axes.set_ylabel(r'10 m wind speed [m/s]',fontsize=FS,color=LC)
axes.plot(time,speed,color=LC,linewidth=1.5)
axes.plot([time[0],time[-1]],[mean,mean],color=LC,linewidth=0.5,linestyle='--')
axis_lims = [time[0], 1.0*time[-1], 0., 20.0]
axes.axis(axis_lims)
#title_string = r"p3_c2    6.25 iper nested run:     $<u_{10}>=12$ m/s      $<\theta>=\pi/4$ "
#axes.set_title(title_string,fontsize=FS)
axes.tick_params(axis='y', labelcolor=LC)
axes.set_yticks([0,5,10,15,20])


mean = 1./4.   # mean diriction/pi
LC='b'
axes2 = axes.twinx()
axes2.tick_params(direction='out', top=False, right=True,) # Turn ticks out
axes2.tick_params(axis='both', which='major', labelsize=FS)
axes2.spines['top'].set_visible(False)                     # Get rid of top axis line
axes2.set_ylabel(r'10 m wind direction',fontsize=FS,color=LC)
axes2.plot(time,dir/np.pi,color=LC,linewidth=1.5)
axes2.plot([time[0],time[-1]],[mean,mean],color=LC,linewidth=0.5,linestyle='--')
#axes2.plot(time,0.25*np.ones_like(time),'--',color=LC,linewidth=0.5)
axis_lims = [time[0], 1.0*time[-1], -0.5, 1]
axes2.axis(axis_lims)
axes2.tick_params(axis='y', labelcolor=LC)

axes2.set_yticks([-.5,0,.5,1])
axes2.set_yticklabels(['$-\pi/2$','0','$\pi/2$','$\pi$'])

plt.savefig(plot_file,dpi=300)      # save plot file







ax = plt.subplot(111, projection='polar')
inc=60*1    # 1 hrs for dt=1 min
ax.plot(dir[0:-1:inc], speed[0:-1:inc],'k.',markersize=0.25)
ax.set_rmax(20)
ax.set_rticks([5, 10, 15, 20])  # less radial ticks
ax.set_rlabel_position(225)  # get radial labels away from plotted line
ax.grid(True)
plt.show()
#plt.savefig(plot_file,dpi=300)      # save plot file
