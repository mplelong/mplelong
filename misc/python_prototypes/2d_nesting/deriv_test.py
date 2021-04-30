#!/usr/local/bin/python
"""
Test driver for testing the differentiation scheme based 
Bernoulli/cosine approach
"""

#---------------------------------------------------------------------------------------
#  import and rename the various modules to be used
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from solver_utilities import evaluate_basis_functions,prepare_BernoulliCosine,ddx 
                                                  
#-----------------------------------------------------------------------------------------


#---------------------------------------------------------------------------
#  set the sample problem
#---------------------------------------------------------------------------
L = np.pi
nx=257
x = np.linspace(0,L,nx)

K=2.; f = np.sin(K*x) ; f_x = K*np.cos(K*x);
title_string = r" exact (r) and computed (b) derivatives      $f=\sin 2x$"
title_string2 = r" log10 error in derivative estimate    $f=\sin 2x$"
 
alpha=L/8 ; f = np.exp(-alpha*x); f_x = -alpha*f 
title_string = r" exact (r) and computed (b) derivatives      $f=e^{-\alpha x}~,~ \alpha=L/8$"
title_string2 = r" log10 error in derivative estimate    $f=e^{-\alpha x}~,~ \alpha=L/8$"


xnorm1 = np.max(np.abs(f_x))


#----------------------------------------------------
#  setup the Bernoulli/Cosine differentiation scheme
#----------------------------------------------------
frac = 0.075
Q = 9            # 5 terms in the singular series

# build the basis functions and their derivs expanded about the 2 endpoints
basis_functions = evaluate_basis_functions(x,Q)

# build and factor expansion coefficient matrices for Bernoulli/Cosine scheme
LU = prepare_BernoulliCosine(x,Q)

#-----------------------------------------------------------
# compute the derivative
#-----------------------------------------------------------
df = ddx(f,x,frac,Q,LU,basis_functions)



#---------------------------------------------------------------------------
#  quantify the errors
#---------------------------------------------------------------------------
error1 = np.log10(np.abs(f_x-df)/xnorm1)


#---------------------------------------------------------------------------------------
#  set up the figure for plotting exact and computed derivatives
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r"$f'$",fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False) 
axes.set_title(title_string,fontsize=12)

axes.plot(x/L,f_x,'r',x/L,df,'bo',markersize=1)
axes.grid()
		
		
plt.savefig('deriv_test.eps',dpi=300,bb_inches='tight')      # save plot file


#---------------------------------------------------------------------------------------
#  set up the figure for plotting errors
#---------------------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
axes.set_ylabel(r'',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(True)                     	# Get rid of right axis line
axes.spines['right'].set_visible(False)    
axes.set_title(title_string2,fontsize=12)
axes.plot(x/L,error1)
axes.grid()
		
		
plt.savefig('error.eps',dpi=300,bb_inches='tight')      # save plot file
