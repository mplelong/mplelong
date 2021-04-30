#-----------------------------------------------------------------------------------------
#  Illustrate combined Bernoulli polynomial cosine expansion approach to differentiating
#  discrete data given at equally spaced grid points in the closed interval [0,L]
#    
#    -----------------------------
#      uses the set of functions defined in 
#         Bernoulli_polynomials.py
#         Fourier_stuff.py
#    -----------------------------
#-----------------------------------------------------------------------------------------
import os,math,sys          
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from Bernoulli_polynomials \
 import setup_factor_matrix,solve_for_expansion_coeffs,U_series_expansion,even_extend
    
from Fourier_stuff         \
 import dst,dst_filtered,compute_spectrum,fourier_filter






#-------------------------------------------------------
#  define a discrete domain [0,l] w/ L=1 for simplicity
#-------------------------------------------------------
nx=256 ;  L=1.0 ; dx = L/(nx-1.); P=2.*L
x = np.linspace(0.,L,nx); 
inc=3   # inc for showing exact deriv with visible markersize



#---------------------------------------------------
#    define a discrete test function on [0,L]=[0,1]
#---------------------------------------------------
gamma = L/100. ; g = -1*np.exp(-((x-.425*L)/gamma)**2) ; g_x = -(2.*(x-.425*L)/gamma**2) * g
alpha = 1.5 ; f = np.exp(alpha*x) + g ; f_x = alpha*np.exp(alpha*x) + g_x
#K = 2.*np.pi ; phase = np.pi/4. ; f = np.sin(K*x+phase) ; f_x = K*np.cos(K*x+phase)


#------------------------------------------------------------------
#  use standard cosine transform to differentiate f
#  discontinuous derivative should render this estimate "terrible"
#------------------------------------------------------------------
n=1 ; flag=-1
f_x_cos = dst(f,n,L,flag)



#-----------------------------------------------------------------
#  set up matrix problems for coeffs for the 2 series expansions
#-----------------------------------------------------------------
Q = 7      # series containing U_1, U_3, U_5 and U_7
           # i.e. B_2, B_4, B_6, B_8
x0 = 0.                                   # even extension of f has discontinuous deriv at x=0
lu_A,piv_A = setup_factor_matrix(x,Q,x0)
x0 = L                                    # even extension of f also has discontinuous deriv at x=L
lu_B,piv_B = setup_factor_matrix(x,Q,x0)

#---------------------------------------------------------------
#  solve the matrix problems for coeffs for 2 series expansions
#---------------------------------------------------------------
x0 = 0.
A = solve_for_expansion_coeffs(lu_A,piv_A,f,Q,x0)
x0 = L
B = solve_for_expansion_coeffs(lu_B,piv_B,f,Q,x0)



#---------------------------------------------------------------
#  construct the 2 series for x in [0,L], add them to get f_s(x)
#  also differentiate this series
#---------------------------------------------------------------

x0 = 0. ; coeffs = A
s1,s1_x = U_series_expansion(x,Q,x0,coeffs)

x0 = L ; coeffs = B
s2,s2_x = U_series_expansion(x,Q,x0,coeffs)

# combine the two series
f_s = s1 + s2 
f_s_x = s1_x + s2_x  

# extract the smooth, even extendable part that is well approximated by a cosine series
f_Q = f - f_s    # should be Q times differentiable when even expanded 




#---------------------------------------------------------------
#  use standard cosine transform to differentiate f_Q
#---------------------------------------------------------------
n=1 ; flag=-1
#f_Q_x = dst(f_Q,n,L,flag)
f_Q_x = dst_filtered(f_Q,n,L,flag,0.05)
f_x_computed = f_Q_x + f_s_x
error = np.log10( np.abs(f_x - f_x_computed)/np.max(np.abs(f_x)) )

#----------------------------------------------------------
# define the [0,2L) extended domain and construct
# the even extension of f and call it F
# F has periodicity P=2L 
#----------------------------------------------------------
M = (nx-1)*2
x_ext = np.linspace(0,2*L-dx,M)
F = even_extend(f)
F_Q = even_extend(f_Q)
F_S = even_extend(f_s)



#---------------------------------------------------------------
#  compute wavenumber spectrum of f and f_Q
#---------------------------------------------------------------
flag = -1    # Fourier analyze even extension of f
k,S = compute_spectrum(f,flag,L)
k,S_Q = compute_spectrum(f_Q,flag,L)

#  frac < 0 ==> 2/3 rule filter,  frac=0 for no filtering
frac = 0.; 
filter = fourier_filter(k,frac) ; kf = k*filter


#-------------------------------------------------------------------------------------------
#  set up the figure for plotting spectra of f(x) and f_Q(x)
#-------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(6,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.15, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$k$',fontsize=14)
axes.set_ylabel(r'$S(k)$',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     	# Get rid of top axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"$S(k)~~~{\rm and}~~~S_Q(k)~~~$  nx=%d $~~~$  Q=%d"  %(nx,Q)   
axes.set_title(title_string,fontsize=14)


#---------------------------------------------------------------------------------------
#   plot the data 
#---------------------------------------------------------------------------------------
inc=1
plt.loglog(k[1:-1:inc],S[1:-1:inc]*filter[1:-1:inc],'k',markersize=2)
#plt.loglog(k[1:-1:2],S[1:-1:2]*kf[1:-1:2],'k+')
plt.loglog(k,S_Q*filter,'r')
#plt.loglog(k,S_Q*kf,'r:')
#plt.loglog(k,S_Q*kf**6,'r')



#---------------------------------------------------------------------------------------
#   impose axis limits, add tick marks 
#---------------------------------------------------------------------------------------
#axis_lims = [0.,1., min, max]  
#plt.axis(axis_lims)

# x axis   
#majorLocator = MultipleLocator(0.25)
#majorFormatter = FormatStrFormatter('%.2f')    # %d
#minorLocator = MultipleLocator(.05)
#axes.xaxis.set_major_locator(majorLocator)
#axes.xaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
#axes.xaxis.set_minor_locator(minorLocator)

# y axis    
#majorLocator = MultipleLocator(1)
#majorFormatter = FormatStrFormatter('%d')
#minorLocator = MultipleLocator(.25)
#axes.yaxis.set_major_locator(majorLocator)
#axes.yaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
#axes.yaxis.set_minor_locator(minorLocator)


#---------------------------------------------------------------------------------------
#   save the figure 
#--------------------------------------------------------------------------------------- 
fname = 'spectra.eps' 
plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution








#-------------------------------------------------------------------------------------------
#  set up the figure for plotting the test function and the 2 series expansions separately
#-------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(6,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
#axes.set_ylabel(r'$u$',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     	# Get rid of top axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"$f(x),~s_0(x),~s_L(x)~~~$  nx=%d $~~~$  Q=%d"  %(nx,Q)   
axes.set_title(title_string,fontsize=14)


#---------------------------------------------------------------------------------------
#   plot the data 
#---------------------------------------------------------------------------------------
plt.plot(x/L,f,'k',linewidth=2) 
plt.plot(x/L,s1,'b--',x,s2,'r--',linewidth=1)
#plt.grid()



#---------------------------------------------------------------------------------------
#   impose axis limits, add tick marks 
#---------------------------------------------------------------------------------------
max = 1.2*np.max( [np.max(np.abs(s1)), np.max(np.abs(s2))] ) ;  min = -max 
axis_lims = [0.,1., min, max]  
plt.axis(axis_lims)

# x axis   
majorLocator = MultipleLocator(0.25)
majorFormatter = FormatStrFormatter('%.2f')    # %d
minorLocator = MultipleLocator(.05)
axes.xaxis.set_major_locator(majorLocator)
axes.xaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
axes.xaxis.set_minor_locator(minorLocator)

# y axis    
majorLocator = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(.25)
axes.yaxis.set_major_locator(majorLocator)
axes.yaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
axes.yaxis.set_minor_locator(minorLocator)


#---------------------------------------------------------------------------------------
#   save the figure 
#--------------------------------------------------------------------------------------- 
fname = 'f_s0_sL_vs_x.eps' 
plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution




#-------------------------------------------------------------------------------------------
#  set up the figure for plotting the test function f, f_s and f_Q, all even-extended
#-------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(6,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
axes = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
axes.set_xlabel(r'$x/L$',fontsize=14)
#axes.set_ylabel(r'$u$',fontsize=14)
axes.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
axes.tick_params(axis='both', which='major', labelsize=12)
axes.spines['top'].set_visible(False)                     	# Get rid of top axis line
axes.spines['right'].set_visible(False) 
    
title_string = r"$f(x),~f_s(x),~f_Q(x)~~~$  nx=%d $~~~$  Q=%d"  %(nx,Q)   
axes.set_title(title_string,fontsize=14)


#---------------------------------------------------------------------------------------
#   plot the data 
#---------------------------------------------------------------------------------------
plt.plot(x_ext/L,F,'k',linewidth=2) 
plt.plot(x_ext/L,F_Q,'r',linewidth=2)
plt.plot(x_ext/L,F_S,'k--',linewidth=1)
#plt.grid()



#---------------------------------------------------------------------------------------
#   impose axis limits, add tick marks 
#---------------------------------------------------------------------------------------
max = 1.2*np.max( [np.max(np.abs(f)), np.max(np.abs(f_s)), np.max(np.abs(f_Q)) ] ) 
min = -max 
axis_lims = [0.,2., min, max]  
plt.axis(axis_lims)

# x axis   
majorLocator = MultipleLocator(0.25)
majorFormatter = FormatStrFormatter('%.2f')    # %d
minorLocator = MultipleLocator(.05)
axes.xaxis.set_major_locator(majorLocator)
axes.xaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
axes.xaxis.set_minor_locator(minorLocator)

# y axis    
majorLocator = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(.25)
axes.yaxis.set_major_locator(majorLocator)
axes.yaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
axes.yaxis.set_minor_locator(minorLocator)


#---------------------------------------------------------------------------------------
#   save the figure 
#--------------------------------------------------------------------------------------- 
fname = 'f_f_s_f_Q_vs_x.eps' 
plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution




#-------------------------------------------------------------------------------------------
#  set up the figure for plotting deriv and its estimate and the error
#-------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(6,4),dpi=300)              		# fig size in inches
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
ax1 = fig.add_axes([0.10, 0.15, 0.75, 0.75])    			# lower, bottom, width, height in ( 0 to 1)
ax1.set_xlabel(r'$x/L$',fontsize=14)
#axes.set_ylabel(r'$u$',fontsize=14)
ax1.tick_params(direction='out', top=False, right=False,) 	# Turn ticks out
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.spines['top'].set_visible(False)                     	# Get rid of top axis line
ax1.spines['right'].set_visible(False) 
    
title_string = r"$f'(x),~{\rm computed~derivative}~~~$  nx=%d $~~~$  Q=%d"  %(nx,Q)   
ax1.set_title(title_string,fontsize=14)





#---------------------------------------------------------------------------------------
#   plot the data 
#---------------------------------------------------------------------------------------
ax1.plot(x/L,f_x_computed,'r',linewidth=2)
ax1.plot(x/L,f_x_cos,'b',linewidth=.5)
ax1.plot(x[0:-1:inc]/L,f_x[0:-1:inc],'ko',markersize=3)

ax2=ax1.twinx()
ax2.tick_params(direction='out', top=False, right=True,) 	# Turn ticks out
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.spines['top'].set_visible(False)                     	# Get rid of top axis line
ax2.spines['right'].set_visible(True) 

ax2.plot(x/L,error,'k',x/L,error,'k.',linewidth=0.5,markersize=2)
#ax2.set_ylabel("log10 error",color="k",fontsize=14)

#---------------------------------------------------------------------------------------
#   impose axis limits, add tick marks 
#---------------------------------------------------------------------------------------
if( np.min(f_x) < 0 ): min = 1.5*np.min(f_x) ; 
if( np.min(f_x) > 0 ): min = .75*np.min(f_x) ;
max = 1.5*np.max(f_x) 
axis_lims = [0.,1.0, min, max]  
ax1.axis(axis_lims)

# x axis   
majorLocator = MultipleLocator(0.25)
majorFormatter = FormatStrFormatter('%.2f')    # %d
minorLocator = MultipleLocator(.05)
ax1.xaxis.set_major_locator(majorLocator)
ax1.xaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
ax1.xaxis.set_minor_locator(minorLocator)

# y axis  for ax1  
majorLocator = MultipleLocator(25)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(5)
ax1.yaxis.set_major_locator(majorLocator)
ax1.yaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
ax1.yaxis.set_minor_locator(minorLocator)


min = -14 ; max = 0 
axis_lims = [0.,1.0, min, max]  
ax2.axis(axis_lims)

# y axis  for ax2  
majorLocator = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(.25)
ax2.yaxis.set_major_locator(majorLocator)
ax2.yaxis.set_major_formatter(majorFormatter)
# for the minor ticks, use no labels; default NullFormatter
ax2.yaxis.set_minor_locator(minorLocator)



#---------------------------------------------------------------------------------------
#   save the figure 
#--------------------------------------------------------------------------------------- 
fname = 'f_derivs.eps' 
plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution

