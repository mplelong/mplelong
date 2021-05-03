#!python
"""
Routine to read in 2d YZ slice data, distributed or concatenated,
and make an image of normalized v (color) with overlying isopycnals.
"""
#---------------------------------------------------------------------------------------
#  choose one
#---------------------------------------------------------------------------------------
do_print_eps = 0              #  whether to print eps files for each frame
do_print_png = 1              #  whether to print png files for each frame


do_second_scalar = False
do_concat = False              # true if data files have not yet been concatenated

#---------------------------------------------------------------------------------------
#  import and name the various modules I'll use
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os, math and system modules
import numpy as np
import matplotlib
#  THIS BACKEND WILL AVOID ANY CONNECTIONS TO A REMOTE X SERVER
matplotlib.use('Agg') # Must be after importing matplotlib and before importing matplotlib.pyplot or pylab!
#
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.io import netcdf
#-----------------------------------------------------------------------------------------
plt.close('all')


#-----------------------------------------------------------------------------------------
#  run specific parameters
#-----------------------------------------------------------------------------------------
g = 9.81                        # [m/s2]
lat = 53.43906*2.0*np.pi/360.0  # [radians]    child 2 grid
Omega = 2.0*np.pi/(24.*3600.)   # [rad/sec]
f0 = 2.*Omega*np.sin(lat)       # [1/s]   ~55N
iper = 2*np.pi/f0               # [s]
Lz = 3500.                      # [m]
rho0 = 1027.                    # [kg/m3]    density as depth --> infty

do_second_scalar = False

slice_dir = '/../../output/slices/2D/'
data_dir = '../../output/2D/'
plot_dir = '/../../output/figures/YZ_1/'
fn_root = 'YZ_1_'                     # YZ_0, YZ_1 and YZ_2 are the 3 planes saved

k_indices = [540,530,510,499,488,453] # k indices of XY planes saved, can be marked on figure
add_contours = False

p1 = 8                                #  number of processors splitting ny and nz
p2 = 20                               #  number of processors splitting ny and nz
iproc = -9999                         # 0 <= iproc < p1  only used if do_concat is true

#-----------------------------------------------------------------------------------------
#  time slices to be processed
#-----------------------------------------------------------------------------------------
slices = np.arange(0,1601,1)        #  list of all time slices to process   
last_slice=slices[-1]

if not os.path.exists(slice_dir):
    cmd = 'mkdir -p ' + slice_dir
    os.system(cmd)

if not os.path.exists(plot_dir):
    cmd = 'mkdir -p ' + plot_dir
    os.system(cmd)


frame=0        
for islice in slices:
    print "processing time slice: %d " % (islice)

#---------------------------------------------------------------------------------------
#  extract slice and concatenate files to create a YZplane snapshot  NB iproc=0
#---------------------------------------------------------------------------------------
    if( do_concat ):
     cmd = "./concat_YZ.pl " + str(islice) + " " + str(iproc) + " " + str(p2) + " " + data_dir + " YZ " + slice_dir +  " YZ  0"
     print cmd
     os.system(cmd)

#---------------------------------------------------------------------------------------
#  read in data from the desired netcdf file
#---------------------------------------------------------------------------------------
    fname = slice_dir + fn_root + str(islice).zfill(6) + '.nc'
    ncf = netcdf.netcdf_file(fname, 'r')
    y = ncf.variables['y'].data.copy()         # [m]
    z = ncf.variables['z'].data.copy()         # [m]
    u = ncf.variables['u'].data.copy()         # [m/s]
    v = ncf.variables['v'].data.copy()         # [m/s]
    w = ncf.variables['w'].data.copy()         # [m/s]
    rho = ncf.variables['s1'].data.copy()      # [kg/m3]
    if( do_second_scalar ):
        dye = ncf.variables['s2'].data.copy()  # [1]
    t = ncf.variables['time'].data.copy()      # [s]
    ncf.close()
    Lz = z[-1]-z[0]                         # [m]
    Ly = y[-1]-y[0]                         # [m]
    Lx = 1.0                                # [m]
    nx = 1
    ny = y.size
    nz = z.size
    dz = z[1]-z[0] ; dy=y[1]-y[0]          # [m]
#---------------------------------------------------------------------------------------
    
    u = u.reshape(nz,ny,nx).squeeze()
    v = v.reshape(nz,ny,nx).squeeze()
    w = w.reshape(nz,ny,nx).squeeze()
    rho = rho.reshape(nz,ny,nx).squeeze()
    if( do_second_scalar ):
        dye = dye.reshape(nz,ny,nx).squeeze()
    ndims = np.size( v.shape )
    Y,Z = np.meshgrid(y,z)
    
    

#---------------------------------------------------------------------------------------
#  generic setup for the figure
#---------------------------------------------------------------------------------------
    FS = 10                                                      # font size
    DPI=600
    fig = plt.figure(figsize=(5,5),dpi=DPI)                 # fig size in inches
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    axes = fig.add_axes([0.15, 0.15, 0.8, 0.8])               # 
    axes.tick_params(direction='out', top=True, right=True,) # Turn ticks out
    axes.tick_params(axis='both', which='major', labelsize=FS)
    axes.spines['top'].set_visible(True)                        # Get rid of or keep top axis line
    axes.spines['right'].set_visible(False)

#---------------------------------------------------------------------------------------
# define domain extent for image map, axis limits for plot and tick marks
#---------------------------------------------------------------------------------------
    axis_lims = [0., Ly/1000., 0., Lz]
    extent = [y[0]/1000., y[-1]/1000., z[0], z[-1]]

    XmajorLocator   = MultipleLocator(10)            # major ticks x
    XmajorFormatter = FormatStrFormatter('%d')       # %d is integer  %.2f
    XminorLocator   = MultipleLocator(2)             # minor ticks x

    
    YmajorLocator   = MultipleLocator(250)
    YmajorFormatter = FormatStrFormatter('%d')
    YminorLocator   = MultipleLocator(50)
    yticks = np.arange(0,3500+1,250)
    axes.set_yticks(yticks)
    axes.set_yticklabels(yticks[::-1])

    w_min=-1.5; w_max=1.5
    CMAP = plt.get_cmap('gray')    # bwr   gray
    
    extent = [y[0]/1000., y[-1]/1000., z[0], z[-1]]   
    im = plt.imshow(w*1000,extent=extent,interpolation='bilinear',origin='lower',cmap=CMAP,aspect=0.0625/2.)
    im.set_clim(w_min,w_max)
    
    cb = plt.colorbar(orientation='vertical', shrink=0.9, pad = 0.05) 
    cb.set_ticks( np.linspace(w_min,w_max,5) )
    cb.set_ticklabels( np.linspace(w_min,w_max,5) )
    cb.ax.tick_params(axis='y', which='major', labelsize=FS-2)
    cb.ax.set_title('[mm/s]',fontsize=FS-1)
    
    
    
#---------------------------------------------------------------------------------------
#  color contour plot for rho
#   CMAP = plt.get_cmap('hot')    # bwr
#---------------------------------------------------------------------------------------
    if(add_contours):
    	ncontours = 64 ; levels = np.linspace(1021,1026,ncontours)
    	plt.contour(Y/1000.,Z,rho,levels,colors='r',extend='both',extent=extent,zorder=2,linewidths=0.25)   # '0.50'

#---------------------------------------------------------------------------------------
# add markers showing depths of saved XY planes
#---------------------------------------------------------------------------------------
    yy=np.asarray([y[-1],y[-1],y[-1],y[-1],y[-1],y[-1]])
    zz=np.asarray([z[k_indices[0]],z[k_indices[1]],z[k_indices[2]],z[k_indices[3]],z[k_indices[4]],z[k_indices[5]]])
    plt.plot(yy/1000.,zz,"<",markersize=4,markerfacecolor="k",markeredgecolor="k") 
    
    yy=np.asarray([y[0],y[0],y[0],y[0],y[0],y[0]])
    plt.plot(yy/1000.,zz,">",markersize=4,markerfacecolor="k",markeredgecolor="k")           
    
    
#---------------------------------------------------------------------------------------
# add tick marks, impose axis limits
#---------------------------------------------------------------------------------------
    plt.axis(axis_lims)

    axes.xaxis.set_major_locator(XmajorLocator)
    axes.xaxis.set_major_formatter(XmajorFormatter)
    #for the minor ticks, use no labels; default NullFormatter
    axes.xaxis.set_minor_locator(XminorLocator)

    #axes.yaxis.set_major_locator(YmajorLocator)
    #axes.yaxis.set_major_formatter(YmajorFormatter)
    #for the minor ticks, use no labels; default NullFormatter
    axes.yaxis.set_minor_locator(YminorLocator)
    
    #plt.gca().invert_yaxis()

#---------------------------------------------------------------------------------------
#  axis labels and title
#---------------------------------------------------------------------------------------
    axes.set_xlabel(r'$y \, \rm{[km]}$',fontsize=FS+2)
    axes.set_ylabel(r'depth \, \rm{[m]}',fontsize=FS+2)
    title_string = r" t=%.2f inertial periods " %(t/iper)
    axes.set_title(title_string, fontsize=FS+2)

#---------------------------------------------------------------------------------------
#  possibly, print the figure as an eps file
#---------------------------------------------------------------------------------------
    if do_print_eps == 1:
        fname = plot_dir + fn_root + 'w_' + str(frame).zfill(4) + '.eps'
        plt.savefig(fname,dpi=DPI,bb_inches='tight') # save an eps w/ nice resolution

#---------------------------------------------------------------------------------------
#  possibly, print the figure as a png file
#---------------------------------------------------------------------------------------
    if do_print_png == 1:
        fname = plot_dir + fn_root + 'w_' + str(frame).zfill(4) + '.png'
        plt.savefig(fname,dpi=DPI,bb_inches='tight') # save a png

    frame = frame+1
    plt.close(fig)
    del axes
    plt.close('all')
    
exit()


    
