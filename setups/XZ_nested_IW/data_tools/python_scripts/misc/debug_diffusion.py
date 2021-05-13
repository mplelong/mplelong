"""
Script to look at diffusion debug data
"""
import os,math,sys              
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

data_dir = '/Users/kraig/flow_solve_BC/output/debug_data/'
Q = 9  ; M=(Q+1)/2 ; filter_fraction=0.00
Lx =100.e3 ; Ly = 30.e3 ; Lz = 600

#----------------------------------------------------------------------
#   look at x diffusion data
# write(1,*) ig,kx(ig),kxfilter(ig),diff_factor(j,ig,1),tmpY(j,ig,k,1)
#----------------------------------------------------------------------
fn = data_dir + 'x_diffusion'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
k  = df[:,1]
kf = df[:,2]
diff_factor =df[:,3]
u_hat =df[:,4]
n = k.size
dk = k[1]
L = np.pi/dk

plt.subplot(2,1,1)
plt.plot( k,np.log10(np.abs(u_hat)),'b' )
plt.plot( k,np.log10(np.abs(u_hat)*diff_factor),'r',linewidth=0.5 )
plt.xlabel('x wavenumber k')
plt.title(' abs u_hat(k)')

plt.subplot(2,1,2)
plt.plot( k,diff_factor,'b' )
plt.xlabel('x wavenumber k')
plt.title('diffusion factor')

plotfile = data_dir + 'x_diffusion.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300)
plt.close()


#----------------------------------------------------------------------
#   look at y diffusion data
# j,ky(j),kyfilter(j),diff_factor(j,i,1),tmpY(j,i,k,1)
#----------------------------------------------------------------------
fn = data_dir + 'y_diffusion'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
k  = df[:,1]
kf = df[:,2]
diff_factor =df[:,3]
u_hat =df[:,4]
n = k.size
dk = k[1]
L = np.pi/dk

plt.subplot(2,1,1)
plt.plot( k,np.log10(np.abs(u_hat)),'b' )
plt.plot( k,np.log10(np.abs(u_hat)*diff_factor),'r',linewidth=0.5 )
plt.xlabel('y wavenumber k')
plt.title(' abs u_hat(k)')

plt.subplot(2,1,2)
plt.plot( k,diff_factor,'b' )
plt.xlabel('y wavenumber k')
plt.title('diffusion factor')

plotfile = data_dir + 'y_diffusion.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300)
plt.close()


#----------------------------------------------------------------------
#   look at z diffusion data
# k,kz(k),kzfilter(k),diff_factor(k,1),wrk(k)
#----------------------------------------------------------------------
fn = data_dir + 'z_diffusion'   
df = pd.read_csv(fn,delim_whitespace=True,header=None).to_numpy()
k  = df[:,1]
kf = df[:,2]
diff_factor =df[:,3]
u_hat =df[:,4]
n = k.size
dk = k[1]
L = np.pi/dk

plt.subplot(2,1,1)
plt.plot( k,np.log10(np.abs(u_hat)),'b' )
plt.plot( k,np.log10(np.abs(u_hat)*diff_factor),'r',linewidth=0.5 )
plt.xlabel('z wavenumber k')
plt.title(' abs u_hat(k)')

plt.subplot(2,1,2)
plt.plot( k,diff_factor,'b' )
plt.xlabel('z wavenumber k')
plt.title('diffusion factor')

plotfile = data_dir + 'z_diffusion.eps'
print(' ....... saving file ', plotfile )
plt.savefig(plotfile,dpi=300)
plt.close()




