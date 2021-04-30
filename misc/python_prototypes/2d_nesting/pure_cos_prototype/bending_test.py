#!/usr/local/bin/python
"""
Test driver for testing the "diffusion" approach to bending functions near the
boundary enough so that their derivatives can be calculated accurately in the interior
using term by term differentiation of the cos transform. The idea is then that, the
interior solution can then be extrapolated back to the boundary giving a "usable"
global approximation to the true derivative.
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from solver_utilities import diffuse_near_boundary,ddx,fd_deriv
                                                  
#-----------------------------------------------------------------------------------------


#---------------------------------------------------------------------------
#  set the sample problem
#---------------------------------------------------------------------------
L = np.pi
nx = 257
x = np.linspace(0,L,nx)
dx = x[1]-x[0]
npts = 9         # width parameter for boundary region
frac = 0.075     # filter fraction for cos differentiation


#----------------------------------------------------
#  3 test functions and derivatives to try
#----------------------------------------------------
K=8.; f = np.sin(K*x) ; exact_deriv = K*np.cos(K*x) ; xx=1.2; yy=1.2    #(factors for plot limits)
#alpha=3./L ; f=np.exp(alpha*x); exact_deriv = alpha*f ; xx=0.5; yy=1.2
#alpha=-0.5/L ; f=np.exp(alpha*x); exact_deriv = alpha*f ; xx=1.2; yy=0.8

xnorm = np.max(np.abs(exact_deriv))

#-------------------------------------------------------------------------------------
# (1) use finite differences to compute the deriv at a few points near the boundaries
#-------------------------------------------------------------------------------------
nvals = 3*npts
df_left,df_right = fd_deriv(f,x,nvals)

#----------------------------------------------------
# (2)  bend the given function f near the boundaries
#----------------------------------------------------
fs = diffuse_near_boundary(f,x,npts)

#-----------------------------------------------------------
# (3) compute the derivative with standard cos method
#-----------------------------------------------------------
df = ddx(fs,x,frac)

#-----------------------------------------------------------
# (4) replace end values with finite difference estimates
#-----------------------------------------------------------
df[0:nvals] = df_left[0:nvals]
df[nx-nvals:nx] = df_right[0:nvals]

error = (np.abs(exact_deriv-df)/xnorm)

#---------------------------------------------------------------------------------------
#  create figure and plot original and "bent" function f
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r"$f$",fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False) 
axes.set_title(r"$f(x)~~~{\rm and}~~~f_s(x)$ ",fontsize=12)
axes.plot(x/L,f,'k',x/L,fs,'b',markersize=1)
axes.grid()

axis_lims = [0., 1., np.min(f), np.max(f)]
plt.axis(axis_lims)

plt.savefig('bending_test.eps',dpi=300,bb_inches='tight')      # save plot file


#---------------------------------------------------------------------------------------
#  create figure and plot exact and computed derivatives
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r"$f_x$",fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False) 
axes.set_title(r"$f_x(x)~~~{\rm and}~~~{f_s}_x(x)$ ",fontsize=12)
axes.plot(x/L,exact_deriv,'k',x/L,df,'b',markersize=1)
axes.grid()

axis_lims = [0., 1., xx*np.min(exact_deriv), yy*np.max(exact_deriv)]
plt.axis(axis_lims)

plt.savefig('bending_test_derivs.eps',dpi=300,bb_inches='tight')      # save plot file


#---------------------------------------------------------------------------------------
#  create figure and plot error in derivative estimation
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r"error",fontsize=12)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False) 
axes.set_title(r"$\log_{10}|\epsilon|$",fontsize=14)
axes.semilogy(x/L,error,'k')
axes.grid()

#axis_lims = [0., 1., xx*np.min(exact_deriv), yy*np.max(exact_deriv)]
#plt.axis(axis_lims)

plt.savefig('bending_test_error.eps',dpi=300,bb_inches='tight')      # save plot file
exit()

