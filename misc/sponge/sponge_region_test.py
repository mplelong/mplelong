#!/usr/local/bin/python
"""
Test construction of sponge region
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# set up a rectangular domain with dx NE dz
Lx=5.; Lz=1.
nx=257; nz=129
x=np.linspace(0.,Lx,nx) ; dx = Lx/(nx-1.)
z=np.linspace(0.,Lz,nz) ; dz = Lz/(nz-1.)
f = np.zeros_like(x)
pi = np.pi


	
def tophat(x,x0,x1,beta_0,beta_1):
#---------------------------------------
#   simple smooth tophat function
#---------------------------------------
	import numpy as np
	p = 4
	if( x >= x0 and x <= x1 ):
		T = 1.
	elif( x < x0):
		T = np.exp(-((x-x0)/beta_0)**p)
	elif( x > x1 ):
		T = np.exp(-((x-x1)/beta_1)**p)
	return T
	
def sponge_weight(x,Lx,offset,width):
	import numpy as np
	for i in range(nx):
		x0 = offset*Lx ; x1 = x0 + width*Lx
		beta_0 = (width*Lx)*.05 ; beta_1 = (width*Lx)*.15
		f = tophat(x,x0,x1,beta_0,beta_1)
		
		x1 = (1.-offset)*Lx ; x0 = x1 - width*Lx
		beta_1 = (width*Lx)*.05 ; beta_0 = (width*Lx)*.15
		f = f + tophat(x,x0,x1,beta_0,beta_1)
	return f		
		


offset = .025   # fraction of domain to offset sponge region
width  = 0.05   # fraction of domain for sponge width at each end
for i in range(nx):
	f[i] = sponge_weight(x[i],Lx,offset,width) 	
#plt.plot( x,f )
#plt.show()
#exit()

S = np.zeros( (nx,nz), dtype=float )

for i in range((nx-1)/2+1):
	print(" i = ",i)
	for k in range((nz-1)/2+1):
	
		a = sponge_weight(x[i],Lx,offset,width)
		b = sponge_weight(z[k],Lz,offset,width)
		if( z[k]< offset*Lz or z[k] > (1.-offset)*Lz ): a=0.
		if( x[i]< offset*Lx or x[i] > (1.-offset)*Lx ): b=0. 
		S[i,k] = np.max([a,b])
		
		ii = (nx-1) - i
		S[ii,k] = S[i,k]
		
		kk = (nz-1) - k
		S[i,kk] = S[i,k]
		S[ii,kk] = S[ii,k]



#---------------------------------------------------------------------------------------
#  generic setup for the figure
#---------------------------------------------------------------------------------------
FS = 12                                                   # font size
fig = plt.figure(figsize=(6,4),dpi=300)                 # fig size in inches
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
axis_lims = [0., Lx, 0., Lz]

XmajorLocator   = MultipleLocator(1)            # major ticks x
XmajorFormatter = FormatStrFormatter('%.1f')       # %d is integer  %.2f
XminorLocator   = MultipleLocator(.5)             # minor ticks x

YmajorLocator   = MultipleLocator(1)
YmajorFormatter = FormatStrFormatter('%.1f')
YminorLocator   = MultipleLocator(.5)

v_min=0; v_max=1
CMAP = plt.get_cmap('hot_r')    # bwr     
ncontours = 64 ; levels = np.linspace(v_min,v_max,ncontours)
plt.contourf(x,z,S.transpose(),levels,cmap=CMAP,extend='both',zorder=1,linewidths=0.75) 
plt.clim(v_min,v_max)

fname = 'Sponge.png'
plt.savefig(fname,dpi=300,bb_inches='tight') # save a png
