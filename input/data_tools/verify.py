"""
Script to read in some of the output/debug_data and make some simple plots
to demonstrate/verify key routines in the flow_solve suite
"""
import os,math,sys              
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

data_dir = '../../output/debug_data/'
Q = 5  ; M=(Q+1)/2 ; filter_fraction=0.05
Lx =100.e3 ; Ly = 30.e3 ; Lz = 600

#----------------------------------------------------------------------
#   look at wavenumber arrays
#----------------------------------------------------------------------
fn = data_dir + 'x_wavenumbers'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
# i,kx(i),kxfilter(i)
k  = df[:,1]
kf = df[:,2]
nx = k.size
if( nx > 1 ):
	dk = k[1]
	L = np.pi/dk
	plt.plot( k,kf,'k' )
	plt.xlabel('x wavenumber k')
	plt.title(' x wavenumber filter')
	plotfile = data_dir + 'kx_kxfilter.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()

fn = data_dir + 'y_wavenumbers'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
# i,kx(i),kxfilter(i)
k  = df[:,1]
kf = df[:,2]
ny = k.size
if( ny > 1 ):
	dk = k[1]
	L = np.pi/dk
	plt.plot( k,kf,'k' )
	plt.xlabel('y wavenumber l')
	plt.title(' y wavenumber filter')
	plotfile = data_dir + 'ky_kyfilter.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()

fn = data_dir + 'z_wavenumbers'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
# i,kx(i),kxfilter(i)
k  = df[:,1]
kf = df[:,2]
nz = k.size
if( nz > 1 ):
	dk = k[1]
	L = np.pi/dk
	plt.plot( k,kf,'k' )
	plt.xlabel('z wavenumber k')
	plt.title(' z wavenumber filter')
	plotfile = data_dir + 'kz_kzfilter.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()


#----------------------------------------------------------------------
#   look at y basis functions
#----------------------------------------------------------------------
if( ny > 1 ):
	if( M >= 1 ):
		fn = data_dir + 'y_bases_1'   
		df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
		# n=1:  y(i), U_y(i,n,1), U_y(i,n,2), dU_y(i,n,1), dU_y(i,n,2)
		y    = df[:,0] ; L = y[-1]-y[0] 
		U_0  = df[:,1]
		U_L  = df[:,2]
		plt.plot( y/L,U_0,'r', y/L,U_L,'b' )
		inc = 4
		#plt.plot( y[0:-1:inc]/L,np.flip(U_0[0:-1:inc]),'r*', y[0:-1:inc]/L,np.flip(U_L[0:-1:inc]),'b*' ,markersize=2)

	if( M >= 2) :
		fn = data_dir + 'y_bases_2'   
		df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
		# n=2:  y(i), U_y(i,n,1), U_y(i,n,2), dU_y(i,n,1), dU_y(i,n,2)
		y    = df[:,0] ; L = y[-1]-y[0] 
		U_0  = df[:,1]
		U_L  = df[:,2]
		plt.plot( y/L,U_0,'r--', y/L,U_L,'b--' )
	
	if( M >= 3) :
		fn = data_dir + 'y_bases_3'   
		df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
		# n=3:  y(i), U_y(i,n,1), U_y(i,n,2), dU_y(i,n,1), dU_y(i,n,2)
		y    = df[:,0] ; L = y[-1]-y[0] 
		U_0  = df[:,1]
		U_L  = df[:,2]
		plt.plot( y/L,U_0,'r:', y/L,U_L,'b:' )

	if( M >= 4) :
		fn = data_dir + 'y_bases_4'   
		df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
		# n=4:  y(i), U_y(i,n,1), U_y(i,n,2), dU_y(i,n,1), dU_y(i,n,2)
		y    = df[:,0] ; L = y[-1]-y[0] 
		U_0  = df[:,1]
		U_L  = df[:,2]
		plt.plot( y/L,U_0,'r', y/L,U_L,'b' ,linewidth=0.5)
	
	if( M >= 5) :
		fn = data_dir + 'y_bases_5'   
		df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
		# n=5:  y(i), U_y(i,n,1), U_y(i,n,2), dU_y(i,n,1), dU_y(i,n,2)
		y    = df[:,0] ; L = y[-1]-y[0] 
		U_0  = df[:,1]
		U_L  = df[:,2]
		plt.plot( y/L,U_0,'r', y/L,U_L,'b',linewidth=0.25 )
	
	plt.title('y basis functions expanded about 0 (red) and Ly (blue)')
	plt.xlabel('y/Ly')

	plotfile = data_dir + 'y_basis_functions.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	

#----------------------------------------------------------------------
#   look at z basis functions
#----------------------------------------------------------------------
if( M >= 1 ):
	fn = data_dir + 'z_bases_1'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	# n=1:  z(i), U_z(i,n,1), U_z(i,n,2), dU_z(i,n,1), dU_z(i,n,2)
	z    = df[:,0] ; L = z[-1]-z[0] 
	U_0  = df[:,1]
	U_L  = df[:,2]
	plt.plot( z/L,U_0,'r', z/L,U_L,'b' )
	inc = 4
	#plt.plot( z[0:-1:inc]/L,np.flip(U_0[0:-1:inc]),'r*', z[0:-1:inc]/L,np.flip(U_L[0:-1:inc]),'b*' ,markersize=2)

if( M >= 2) :
	fn = data_dir + 'z_bases_2'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	# n=2:  z(i), U_z(i,n,1), U_z(i,n,2), dU_z(i,n,1), dU_z(i,n,2)
	z    = df[:,0] ; L = z[-1]-z[0] 
	U_0  = df[:,1]
	U_L  = df[:,2]
	plt.plot( z/L,U_0,'r--', z/L,U_L,'b--' )
	
if( M >= 3) :
	fn = data_dir + 'z_bases_3'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	# n=2:  z(i), U_z(i,n,1), U_z(i,n,2), dU_z(i,n,1), dU_z(i,n,2)
	z    = df[:,0] ; L = z[-1]-z[0] 
	U_0  = df[:,1]
	U_L  = df[:,2]
	plt.plot( z/L,U_0,'r:', z/L,U_L,'b:' )

if( M >= 4) :
	fn = data_dir + 'z_bases_4'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	# n=4:  z(i), U_z(i,n,1), U_z(i,n,2), dU_z(i,n,1), dU_z(i,n,2)
	z    = df[:,0] ; L = z[-1]-z[0] 
	U_0  = df[:,1]
	U_L  = df[:,2]
	plt.plot( z/L,U_0,'r', z/L,U_L,'b',linewidth=0.5 )
	
if( M >= 5) :
	fn = data_dir + 'z_bases_5'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	# n=5:  z(i), U_z(i,n,1), U_z(i,n,2), dU_z(i,n,1), dU_z(i,n,2)
	z    = df[:,0] ; L = z[-1]-z[0] 
	U_0  = df[:,1]
	U_L  = df[:,2]
	plt.plot( z/L,U_0,'r', z/L,U_L,'b',linewidth=0.25 )
	
plt.title('z basis functions expanded about 0 (red) and Lz (blue)')
plt.xlabel('z/Lz')

plotfile = data_dir + 'z_basis_functions.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	


#----------------------------------------------------------------------
#   look at function decomposition  y, S(y) S_y(y)
#----------------------------------------------------------------------
fn = data_dir + 'combined_series'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z = df[:,0] 
S   = df[:,1]                   # both left and right series together
S_z = df[:,2]                   # both left and right deriv series together
f   = df[:,3]                   # the original function to be differentiated

plt.plot( z/Lz,f,'k', z/Lz,S,'b', z/Lz,f-S, 'r' )
plt.title('f(z) (k), the series expansion S(z) (b) and smooth part f_Q(z) (r)')
plt.xlabel(r'$z/L_z$')

plotfile = data_dir + 'function_decomposition.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	



#----------------------------------------------------------------------
#   look at derivative decomposition  y, d/dy f_Q,  d/dy S(y), d/dy f
#----------------------------------------------------------------------
fn = data_dir + 'BCderivs'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z = df[:,0] 
f_Q_z   = df[:,1]                   # computed deriv of smooth part
S_z = df[:,2]                       # computed deriv of series part
dfdz   = df[:,3]                    # the BC computed deriv

plt.plot( z/Lz,f_Q_z,'k', z/Lz,S_z,'b', z/Lz,dfdz, 'r' )
plt.title('d/dz f_Q (k)  d/dz S(z) (b)  d/dz f(z) (r) ')
plt.xlabel(r'$z/L_z$')
plt.grid()

plotfile = data_dir + 'derivative_decomposition.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	


#----------------------------------------------------------------------
#   look at N/S boundary treatment, i.e. boundary treatment applied
#   immediately to the initial conditions as a test
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'check_NS_bdry_scheme'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0] ; 
	u  = df[:,1] ; u =u/np.max(np.abs(u))
	v  = df[:,2] ; v =v/np.max(np.abs(v))
	w  = df[:,3] ; w =w/np.max(np.abs(w))
	s1 = df[:,4] ; s1=s1/np.max(np.abs(s1))

	plt.plot(y/Ly,u,'r',y/Ly,v,'g',y/Ly,w,'b',y/Ly,s1,'k')
	plt.title('N/S boundary treatment applied to normalized ICs')
	plt.xlabel('y/Ly')
	plt.grid()

	plotfile = data_dir + 'NS_boundary_scheme.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	


#----------------------------------------------------------------------
#   look at bottom boundary treatment, i.e. boundary treatment applied
#   immediately to the initial conditions as a test
#----------------------------------------------------------------------
fn = data_dir + 'check_B_bdry_scheme'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
offset = 1.e-1
z  = df[:,0]
u  = df[:,1] ; u =u/np.max(np.abs(u))
v  = df[:,2] ; v =v/np.max(np.abs(v)) + offset
w  = df[:,3] ; w =w/np.max(np.abs(w))
s1 = df[:,4] ; s1=s1/np.max(np.abs(s1)) + offset


plt.plot(z/Lz,u,'r',z/Lz,v,'g',z/Lz,w,'b',z/Lz,s1,'k')
plt.title('Bottom boundary treatment applied to normalized ICs')
plt.xlabel(r'$z/L_z$')
plt.grid()

plotfile = data_dir + 'B_boundary_scheme.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	

#----------------------------------------------------------------------
#   look at top boundary treatment, i.e. boundary treatment applied
#   immediately to the initial conditions as a test
#----------------------------------------------------------------------
fn = data_dir + 'check_T_bdry_scheme'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
offset = 1.e-1
z  = df[:,0] ; 
u  = df[:,1] ; u =u/np.max(np.abs(u))
v  = df[:,2] ; v =v/np.max(np.abs(v)) + offset
w  = df[:,3] ; w =w/np.max(np.abs(w))
s1 = df[:,4] ; s1=s1/np.max(np.abs(s1)) + offset

plt.plot(z/Lz,u,'r',z/Lz,v,'g',z/Lz,w,'b',z/Lz,s1,'k')
plt.title('Top boundary treatment applied to normalized ICs')
plt.xlabel('$z/L_z$')
plt.grid()

plotfile = data_dir + 'T_boundary_scheme.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	

#----------------------------------------------------------------------
#   look at psi_v vs y 
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'psi_v_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0] 
	psi_v  = df[:,1] 

	plt.plot(y/Ly,psi_v,'k',linewidth=1)
	plt.title(r'$\psi_v(y)$')
	plt.xlabel(r'$y/Ly$')
	plt.grid()

	plotfile = data_dir + 'psi_v_of_y.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()

#----------------------------------------------------------------------
#   look at psi_u vs x 
#----------------------------------------------------------------------
if( nx > 1 ):
	fn = data_dir + 'psi_u_of_x'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	x  = df[:,0] 
	psi_u  = df[:,1] 

	plt.plot(x/Lx,psi_u,'k',linewidth=1)
	plt.title(r'$\psi_u(x)$')
	plt.xlabel(r'$x/Lx$')
	plt.grid()

	plotfile = data_dir + 'psi_u_of_x.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()


#----------------------------------------------------------------------
#   look at psi_w vs z 
#----------------------------------------------------------------------
fn = data_dir + 'psi_w_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
psi_w  = df[:,1] 

plt.plot(z/Lz,psi_w,'k',linewidth=1)
plt.title(r'$\psi_w(z)$')
plt.xlabel(r'$z/Lz$ ')
plt.grid()

plotfile = data_dir + 'psi_w_of_z_B.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()

fn = data_dir + 'psi_w_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
psi_w  = df[:,1] 

plt.plot(z/Lz,psi_w,'k',linewidth=1)
plt.title(r'$\psi_w(z)$')
plt.xlabel(r'$z/Lz$ ')
plt.grid()

plotfile = data_dir + 'psi_w_of_z_T.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()



#----------------------------------------------------------------------
#   look at ustar, vstar and wstar vs y before augmenting w/ bdry info
#   and ustar_hat etc after adding psi_u, psi_v and psi_w
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'ustar_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0] 
	ustar  = df[:,1] 
	vstar  = df[:,2]
	wstar  = df[:,3] 
	offset = .02*np.max(np.abs(vstar))

	fn = data_dir + 'ustar_hat_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	ustar_hat  = df[:,1] 
	vstar_hat  = df[:,2] 
	wstar_hat  = df[:,3] 

	plt.plot(y/Ly,vstar,'g',y/Ly,vstar_hat+offset,'g--')
	plt.title(r'${v}_*$ and ${\hat {v}}_*$ (--)  (offset)')
	plt.xlabel(r'$y/Ly$')
	plt.grid()

	plotfile = data_dir + 'vstar_vstar_hat_of_y.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	


#----------------------------------------------------------------------
#   look at ustar, vstar and wstar vs x before augmenting w/ bdry info
#   and ustar_hat etc after adding psi_u, psi_v and psi_w
#----------------------------------------------------------------------
if( nx > 1 ):
	fn = data_dir + 'ustar_of_x'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	x  = df[:,0] 
	ustar  = df[:,1] 
	vstar  = df[:,2]
	wstar  = df[:,3] 
	offset = .02*np.max(np.abs(ustar))

	fn = data_dir + 'ustar_hat_of_x'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	ustar_hat  = df[:,1] 
	vstar_hat  = df[:,2] 
	wstar_hat  = df[:,3] 

	plt.plot(x/Lx,ustar,'g',x/Lx,ustar_hat+offset,'g--')
	plt.title(r'${u}_*$ and ${\hat {u}}_*$ (--)  (offset)')
	plt.xlabel(r'$x/Lx$')
	plt.grid()

	plotfile = data_dir + 'ustar_ustar_hat_of_x.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	



#----------------------------------------------------------------------
#   look at ustar, vstar and wstar vs z before augmenting w/ bdry info
#   and ustar_hat etc after adding psi_u, psi_v and psi_w
#    (near BOTTOM)
#----------------------------------------------------------------------
fn = data_dir + 'ustar_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
ustar  = df[:,1]
vstar  = df[:,2] 
wstar  = df[:,3] 

fn = data_dir + 'ustar_hat_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
ustar_hat  = df[:,1] 
vstar_hat  = df[:,2] 
wstar_hat  = df[:,3] 

plt.plot(z/Lz,wstar,'b',z/Lz,wstar_hat,'b:')
plt.title(r'Near bottom:       ${w}_*$ and ${\hat {w}}_*$ (...)')
plt.xlabel(r'$z/Lz$ ')
plt.grid()

plotfile = data_dir + 'wstar_wstar_hat_of_z_B.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()


#----------------------------------------------------------------------
#   look at ustar, vstar and wstar vs z before augmenting w/ bdry info
#   and ustar_hat etc after adding psi_u, psi_v and psi_w
#    (near TOP)
#----------------------------------------------------------------------
fn = data_dir + 'ustar_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
ustar  = df[:,1]
vstar  = df[:,2] 
wstar  = df[:,3] 

fn = data_dir + 'ustar_hat_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
ustar_hat  = df[:,1] 
vstar_hat  = df[:,2] 
wstar_hat  = df[:,3] 

plt.plot(z/Lz,wstar,'b',z/Lz,wstar_hat,'b:')
plt.title(r'Near top:       ${w}_*$ and ${\hat {w}}_*$ (...)')
plt.xlabel(r'$z/Lz$ ')
plt.grid()

plotfile = data_dir + 'wstar_wstar_hat_of_z_T.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()



#----------------------------------------------------------------------
#   look at div(u*) vs y just prior to the pressure solve
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'div_ustar_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0] 
	div= df[:,1] 

	plt.plot(y/Ly,div,'k')
	plt.title(r'$\nabla \cdot {\vec u}_*$')
	plt.xlabel(r'$y/Ly$')
	plt.grid()

	plotfile = data_dir + 'div_ustar_of_y.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	

#----------------------------------------------------------------------
#   look at div(u*) vs z just prior to the pressure solve
#----------------------------------------------------------------------
fn = data_dir + 'div_ustar_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
div= df[:,1] 

plt.plot(z/Lz,div,'k')
plt.title(r'$\nabla \cdot {\vec u}_*$   (bottom)')
plt.xlabel(r'$z/Lz$')
plt.grid()

plotfile = data_dir + 'div_ustar_of_z_B.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	

fn = data_dir + 'div_ustar_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
div= df[:,1] 

plt.plot(z/Lz,div,'k')
plt.title(r'$\nabla \cdot {\vec u}_*$   (top)')
plt.xlabel(r'$z/Lz$')
plt.grid()

plotfile = data_dir + 'div_ustar_of_z_T.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()


#----------------------------------------------------------------------
#   look at phi vs y 
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'phi_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0] 
	phi= df[:,1] 
	rhs= df[:,2]

	plt.plot(y/Ly,phi,'k')
	plt.title(r'$\phi(y)$')
	plt.xlabel(r'$y/Ly$')
	plt.grid()

	plotfile = data_dir + 'phi_of_y.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	

#----------------------------------------------------------------------
#   look at phi vs z 
#----------------------------------------------------------------------
fn = data_dir + 'phi_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
phi= df[:,1] 
rhs= df[:,2]

plt.plot(z/Lz,phi,'k')
plt.title(r'$\phi(z)$    (bottom)')
plt.xlabel(r'$z/L_z$  ')
plt.grid()

plotfile = data_dir + 'phi_of_z_B.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()

fn = data_dir + 'phi_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
phi= df[:,1] 
rhs= df[:,2]

plt.plot(z/Lz,phi,'k')
plt.title(r'$\phi(z)$    (top)')
plt.xlabel(r'$z/L_z$ ')
plt.grid()

plotfile = data_dir + 'phi_of_z_T.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()


#----------------------------------------------------------------------
#   look at grad phi vs y 
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'grad_phi_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0]
	phi_x = df[:,1] 
	phi_y = df[:,2]
	phi_z = df[:,3]

	plt.plot(y/Ly,phi_x,'r',y/Ly,phi_y,'g',y/Ly,phi_z,'b')
	plt.title(r'$\nabla \phi(y)$')
	plt.xlabel(r'$y/Ly$')
	plt.grid()

	plotfile = data_dir + 'grad_phi_of_y.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	

#----------------------------------------------------------------------
#   look at grad phi vs z 
#----------------------------------------------------------------------
fn = data_dir + 'grad_phi_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
phi_x = df[:,1] 
phi_y = df[:,2]
phi_z = df[:,3]

plt.plot(z/Lz,phi_x,'r',z/Lz,phi_y,'g',z/Lz,phi_z,'b')
plt.title(r'$\nabla \phi(z)$   (bottom)')
plt.xlabel(r'$z/L_z$')
plt.grid()

plotfile = data_dir + 'grad_phi_of_z_B.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	

fn = data_dir + 'grad_phi_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
phi_x = df[:,1] 
phi_y = df[:,2]
phi_z = df[:,3]

plt.plot(z/Lz,phi_x,'r',z/Lz,phi_y,'g',z/Lz,phi_z,'b')
plt.title(r'$\nabla \phi(z)$   (top)')
plt.xlabel(r'$z/L_z$  ')
plt.grid()

plotfile = data_dir + 'grad_phi_of_z_T.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	



#----------------------------------------------------------------------
#   look at pressure corrected velocity vs y 
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'u_projected_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0] 
	u = df[:,1] ; u = u/np.max(np.abs(u))
	v = df[:,2] ; v = v/np.max(np.abs(v))
	w = df[:,3] ; w = w/np.max(np.abs(w))

	plt.plot(y/Ly,u,'r',y/Ly,v,'g',y/Ly,w,'b')
	plt.title(r'pressure projected $\vec{u}(y)$')
	plt.xlabel(r'$y/L_y$')
	plt.grid()

	plotfile = data_dir + 'u_projected_of_y.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	

#----------------------------------------------------------------------
#   look at pressure corrected velocity vs z 
#----------------------------------------------------------------------
fn = data_dir + 'u_projected_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
u = df[:,1] ; #u = u/np.max(np.abs(u))
v = df[:,2] ; #v = v/np.max(np.abs(v))
w = df[:,3] ; #w = w/np.max(np.abs(w))

plt.plot(z/Lz,u,'r',z/Lz,v,'g',z/Lz,w,'b')
plt.title(r'pressure projected $\vec{u}(z)$   (bottom)')
plt.xlabel(r'$z/L_z$ ')
plt.grid()

plotfile = data_dir + 'u_projected_of_z_B.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	

fn = data_dir + 'u_projected_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
u = df[:,1] ; #u = u/np.max(np.abs(u))
v = df[:,2] ; #v = v/np.max(np.abs(v))
w = df[:,3] ; #w = w/np.max(np.abs(w))

plt.plot(z/Lz,u,'r',z/Lz,v,'g',z/Lz,w,'b')
plt.title(r'pressure projected $\vec{u}(z)$   (top)')
plt.xlabel(r'$z/L_z$')
plt.grid()

plotfile = data_dir + 'u_projected_of_z_T.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()



#----------------------------------------------------------------------
#   look at final solns after apply_bcs vs y 
#----------------------------------------------------------------------
if( ny > 1 ):
	fn = data_dir + 'final_u_of_y'   
	df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
	y  = df[:,0] 
	u = df[:,1] ; #u = u/np.max(np.abs(u))
	v = df[:,2] ; #v = v/np.max(np.abs(v))
	w = df[:,3] ; #w = w/np.max(np.abs(w))
	s1= df[:,4] ; #s1 = s1/np.max(np.abs(s1))

	plt.plot(y/Ly,u,'r',y/Ly,v,'g',y/Ly,w,'b',y/Ly,s1,'k')
	plt.title(r'solns at end of time step')
	plt.xlabel(r'$y/L_y$')
	plt.grid()

	plotfile = data_dir + 'final_vals_of_y.eps'
	print(' ....... saving file ', plotfile )
	plt.savefig(plotfile,dpi=300,bb_inches='tight')
	plt.close()	

#----------------------------------------------------------------------
#   look at final solns after apply_bcs vs z 
#----------------------------------------------------------------------
fn = data_dir + 'final_u_of_z_B'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
u = df[:,1] ; #u = u/np.max(np.abs(u))
v = df[:,2] ; #v = v/np.max(np.abs(v))
w = df[:,3] ; #w = w/np.max(np.abs(w))
s1= df[:,4] ; #s1 = s1/np.max(np.abs(s1))

plt.plot(z/Lz,u,'r',z/Lz,v,'g',z/Lz,w,'b',z/Lz,s1,'k')
plt.title(r'final solns at end of time step  (bottom)')
plt.xlabel(r'$z/L_z$ ')
plt.grid()

plotfile = data_dir + 'final_vals_of_z_B.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	

fn = data_dir + 'final_u_of_z_T'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
z  = df[:,0] 
u = df[:,1] ; #u = u/np.max(np.abs(u))
v = df[:,2] ; #v = v/np.max(np.abs(v))
w = df[:,3] ; #w = w/np.max(np.abs(w))
s1= df[:,4] ; #s1 = s1/np.max(np.abs(s1))

plt.plot(z/Lz,u,'r',z/Lz,v,'g',z/Lz,w,'b',z/Lz,s1,'k')
plt.title(r'final solns at end of time step  (top)')
plt.xlabel(r'$z/L_z$')
plt.grid()

plotfile = data_dir + 'final_vals_of_z_T.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300,bb_inches='tight')
plt.close()	
















