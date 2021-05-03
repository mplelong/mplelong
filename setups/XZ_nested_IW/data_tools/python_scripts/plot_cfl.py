#!
"""
    simple plot of cfl data
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
#  THIS BACKEND WILL AVOID ANY CONNECTIONS TO A REMOTE X SERVER
# Must be after importing matplotlib and before importing matplotlib.pyplot or pylab!
matplotlib.use('Agg')
#
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#-----------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#  may need to check paths depending on where "input"  is located relative to "flow_solve"
#------------------------------------------------------------------------------------------
data_file = '../../output/cfl.dat'     # t, cfl_x, cfl_y, cfl_z
plot_dir = '../../output/figures/'
plot_file = plot_dir + 'cfl.eps'

if not os.path.exists(plot_dir):
    cmd = 'mkdir -p ' + plot_dir
    os.system(cmd)


#---------------------------------------
# Open file & read data columns...
#---------------------------------------
f = open(data_file, 'r')
f = np.loadtxt(data_file)
time = f[:,0] ; cfl_x = f[:,1]  ; cfl_y = f[:,2] ; cfl_z = f[:,3]
dt = 60.0




fig = plt.figure(figsize=(8,3.75),dpi=300)              # fig size in inches
#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    # lower, bottom, width, height in ( 0 to 1)
axes.tick_params(direction='out', top=False, right=False,) # Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     # Get rid of top axis line
axes.spines['right'].set_visible(False)

time = time/(24.*3600)   # work in days

plt.plot(time,cfl_x,'r',label='cfl x')
plt.plot(time,cfl_y,'g',label='cfl y')
plt.plot(time,cfl_z,'b',label='cfl z')
plt.plot([time[0],time[-1]],[0.2,0.2],'k--')
plt.legend(loc='lower left')

axis_lims = [time[0], 1.05*time[-1], 0., 0.25]
plt.axis(axis_lims)

axes.set_xlabel(r'$t \, \rm[days]$',fontsize=14)
axes.set_ylabel(r'$(v,w)\,dt/(dy,dz)$',fontsize=14)
title_string = r"low resolution parent_run:          dt=%.2f  [s] " %(dt)
axes.set_title(title_string,fontsize=12)


plt.savefig(plot_file,dpi=300,bb_inches='tight')      # save plot file
