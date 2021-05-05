#!/opt/local/bin/python2.7
"""
Routine to read time slices from a spatially concatenated XY netcdf file (2d)
and make a plot of the density and vorticity fields.
usage example on morgon:  
/usr/bin/mpirun -np 2 python plot_XY_slice.py < /dev/null &> slice_log &
"""
#---------------------------------------------------------------------------------------
#  import and name the various python modules I'll use
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
import matplotlib
#  THIS BACKEND WILL AVOID ANY CONNECTIONS TO A REMOTE X SERVER
# Must be after importing matplotlib and before importing matplotlib.pyplot or pylab!
matplotlib.use('Agg') 
#
import matplotlib.pyplot as plt
import colormaps as cmaps
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.io import netcdf
from mpi4py import MPI
from KBW_FUNCTIONS import dst,vort_z
#-----------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------
#  image format
#---------------------------------------------------------------------------------------
do_print_eps = 1              #  whether to print eps files for each frame
do_print_png = 0              #  whether to print png files for each frame
slice_offset = 0              #  offset for filename indices
plot_dir = '../../figures/'




#-----------------------------------------------------------------------------------------
#  run specific parameters
#-----------------------------------------------------------------------------------------
g = 9.81           # [m/s2]
f0 = 1.191e-4      # [1/s]   55N
iper = 2*np.pi/f0  # [s]
Lz = 3500.         # [m]
rho0 = 1027.       # [kg/m3]    density as depth --> infty

do_second_scalar = False

slice_dir = '../../output/slices/2D/'
data_dir = '../../output/2D/'

p1 = 8                                #  number of processors splitting ny and nz
p2 = 24                               #  number of processors splitting ny and nz

#-----------------------------------------------------------------------------------------
#  time slices to be processed
#-----------------------------------------------------------------------------------------
slices = np.arange(0,993,1)  #np.arange(203,602,1)        #  list of all time slices to process   
last_slice=slices[-1]


#-----------------------------------------------------------------------------------------
#   data format 
#   'distributed'   # e.g. 2D/XY_iii-jjj.nc   iii in [0,p1-1] , jjj in [0,p2-1] concatenated time slices
#   'snapshots'     # e.g. slices/2D/XY_xxxxxx.nc , for slice index xxxxxx
#-----------------------------------------------------------------------------------------
data_format = 'snapshots'         # e.g. 2D/XY.nc , concatenated time slices
jproc = 9999                      # not needed for snapshots format



#-----------------------------------------------------------------------------------------
#   initialize MPI 
#-----------------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
numprocs = MPI.COMM_WORLD.Get_size()    # number of processors executing this script
myid = MPI.COMM_WORLD.Get_rank()        # the processor id of this processor
name = MPI.Get_processor_name()         # myid's actual processor name, i.e. amigo.ucsd.edu
if myid==0: 
 print " plot_XY_slice.py running on  %d processors on %s.\n" % (numprocs,name) 



#-----------------------------------------------------------------------------------------
#   divide up the time slices across the numprocs processors
#   each mpi task makes its own list of slices: my_slices
#-----------------------------------------------------------------------------------------
len = np.float( len(slices) )                             # length of global slice list, as float
local_len = np.int(np.ceil(len/numprocs))                 # nominal length of slice list on each processor
i0 = myid*local_len                                       # starting index for this processors slices
i1 = i0 + local_len                                       # (noninclusive) ending index
i1 = min(i1,last_slice+1)                                 # keep last list(s) from overshooting the end
my_slices = slices[i0:i1]                                 # list of slices assigned to processor myid
comm.Barrier()                                            # all start together...

if( myid==0 ):
	if not os.path.exists(plot_dir):
		cmd = 'mkdir -p ' + plot_dir
		os.system(cmd)

plt.close('all')
comm.Barrier()



for islice in my_slices:
    print "processor %d processing time slice: %d " % (myid,islice)
       
#---------------------------------------------------------------------------------------
#  if data is still distributed ==> spatially concatenate each 
#  time slice separately and store in slices/2D 
#        ./concat_XY.pl tslice p1 jproc data_dir data_root out_dir out_root myid
#---------------------------------------------------------------------------------------
    if( data_format=='distributed' ): 
     cmd = "perl ./concat_XY.pl " + str(islice-slice_offset) + " " + str(p1) + " " + str(jproc) + " " + data_dir[0:-1] + " XY " + slice_dir[0:-1] +  " XY  0"
     print cmd
     os.system(cmd)
     fname = slice_dir + 'XY_surf_' + str(islice-slice_offset).zfill(6) + '.nc'

#---------------------------------------------------------------------------------------
#  if already concatenated, read from the file corresponding to the current time 
#---------------------------------------------------------------------------------------
    if( data_format=='snapshots' ): 
     fname = slice_dir + 'XY_surf_' + str(islice-slice_offset).zfill(6) + '.nc'
     
#---------------------------------------------------------------------------------------
#  if serial data, i.e. 1 XY plane with all time slices
#  just extract the desired time slice store in a tmp file 
#---------------------------------------------------------------------------------------     
    if( data_format=='serial' ):
     cmd = 'ncks -dtimedimension,' +  str(islice-slice_offset) + ' ' + data_dir + 'XY.nc tmp.nc'
     os.system(cmd)
     fname = 'tmp.nc'


    print 'Reading file ', fname    
#---------------------------------------------------------------------------------------
#  read in data from the desired netcdf file
#---------------------------------------------------------------------------------------
    ncf = netcdf.netcdf_file(fname, 'r')
    x = ncf.variables['x'].data.copy()                   # [m]
    y = ncf.variables['y'].data.copy()                   # [m]   
    z = ncf.variables['z'].data.copy()                   # [m]
    u = ncf.variables['u'].data.copy().squeeze().transpose()         # [m/s]
    v = ncf.variables['v'].data.copy().squeeze().transpose()         # [m/s]
    w = ncf.variables['w'].data.copy().squeeze().transpose()         # [m/s] 
    rho = ncf.variables['s1'].data.copy().squeeze().transpose()      # [kg/m3]
    if( do_second_scalar ):
        dye = ncf.variables['s2'].data.copy().squeeze().transpose()  # [1]
    t = ncf.variables['time'].data.copy()                # [s]
    ncf.close()
        
    Lx = x[-1]-x[0]                         # [m]
    Ly = y[-1]-y[0]                         # [m]
    nz = 1
    ny = y.size
    nx = x.size
    ndims = np.size( v.shape )
    Y,X = np.meshgrid(y,x)
    print 'shape of 2d data arrays u,v,w,rho: ', v.shape  #  [ny,nx]
    
    if( data_format=='serial' ): os.system('rm -f tmp.nc')
    
#---------------------------------------------------------------------------------------
#  compute buoyancy
#---------------------------------------------------------------------------------------
    b = -g*(rho - rho0)/rho0     # 0 < b < g*delta_rho/rho0 ~ .033
    delta_b = g*4.0/rho0
    
    delta_b = 0.04    # approximate

    
    
#---------------------------------------------------------------------------------------
#  generic setup for the figure
#---------------------------------------------------------------------------------------
    FS = 8    
    fig = plt.figure(figsize=(4,4),dpi=600)              # fig size in inches
    #plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    axes = fig.add_axes([0.15, 0.15, 0.8, 0.8])    # lower, bottom, width, height in ( 0 to 1)
    axes.tick_params(direction='out', top=False, right=False,) # Turn ticks out
    axes.tick_params(axis='both', which='major', labelsize=FS-2)
    axes.spines['top'].set_visible(True)                    
    axes.spines['right'].set_visible(True) 
    
    
#---------------------------------------------------------------------------------------
#  contour the density range w shading, add contour lines 
#  add colorbar and set its face color to match the background color
#---------------------------------------------------------------------------------------
    #CMAP = plt.get_cmap('plasma')    # bwr
    CMAP=cmaps.plasma                 # vidris, magma, inferno, plasma rainbow
    
    extent = [x[0]/1000., x[-1]/1000., y[0]/1000., y[-1]/1000.]   
    im = plt.imshow(b,extent=extent,interpolation='bilinear',origin='lower',cmap=CMAP,aspect=1.0)
    im.set_clim(0.02,delta_b)

    
    cb = plt.colorbar(orientation='vertical', shrink=0.785, pad = 0.007) 
    cb.set_ticks( np.linspace(0.02,delta_b,5) )
    cb.set_ticklabels( np.linspace(0.02,delta_b,5) )
    cb.ax.tick_params(axis='y', which='major', labelsize=FS-2)
    


#---------------------------------------------------------------------------------------
#   impose axis limits, add tick marks 
#--------------------------------------------------------------------------------------- 
    axis_lims = [0.,Lx/1000., 0., Ly/1000.]  
    plt.axis(axis_lims)
    
    majorLocator = MultipleLocator(100)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(25)
    axes.xaxis.set_major_locator(majorLocator)
    axes.xaxis.set_major_formatter(majorFormatter)
    # for the minor ticks, use no labels; default NullFormatter
    axes.xaxis.set_minor_locator(minorLocator)
    
    majorLocator = MultipleLocator(100)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(25)
    axes.yaxis.set_major_locator(majorLocator)
    axes.yaxis.set_major_formatter(majorFormatter)
    # for the minor ticks, use no labels; default NullFormatter
    axes.yaxis.set_minor_locator(minorLocator)

    axes.set_xlabel(r'$x\, \,km$',fontsize=FS)
    axes.set_ylabel(r'$y\, \,km$',fontsize=FS)
    title_string = r"$b$      depth=%.2f m             t=%.2f  inertial periods" %( (Lz-z),t/iper) 
    axes.set_title(title_string, fontsize=FS-2)


#---------------------------------------------------------------------------------------
#  possibly, print the figure as an eps file
#---------------------------------------------------------------------------------------
    if do_print_eps == 1: 
        fname = plot_dir + 'buoy_' + str(islice).zfill(4) + '.eps' 
        plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution  
        
#---------------------------------------------------------------------------------------
#  possibly, print the figure as a png file
#---------------------------------------------------------------------------------------
    if do_print_png == 1: 
        fname = plot_dir + 'buoy_' + str(islice).zfill(4) + '.png' 
        plt.savefig(fname,dpi=600,bb_inches='tight') # save a png  
        
    
    plt.close(fig)
    del axes
    plt.close('all')
    
    

#---------------------------------------------------------------------------------------
#  now compute and plot the vertical vorticity vort_z = u_y - v_x
#---------------------------------------------------------------------------------------
    zeta = vort_z(u,v,Lx,Ly)
    
#---------------------------------------------------------------------------------------
#  generic setup for the vorticity figure
#---------------------------------------------------------------------------------------
    FS = 8    
    fig = plt.figure(figsize=(4,4),dpi=600)              # fig size in inches
    #plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    axes = fig.add_axes([0.15, 0.15, 0.8, 0.8])    # lower, bottom, width, height in ( 0 to 1)
    axes.tick_params(direction='out', top=False, right=False,) # Turn ticks out
    axes.tick_params(axis='both', which='major', labelsize=FS-2)
    axes.spines['top'].set_visible(True)                    
    axes.spines['right'].set_visible(True) 
    
   
#---------------------------------------------------------------------------------------
#  contour the positive and negative vertical vorticity
#---------------------------------------------------------------------------------------
    CMAP = plt.get_cmap('seismic')    # bwr
    vmax = 1.0   # np.amax( np.abs(zeta) )/f0 ; print 'max vorticity/f0 ', vmax
    
    
    extent = [x[0]/1000., x[-1]/1000., y[0]/1000., y[-1]/1000.]   
    im = plt.imshow(zeta/f0,extent=extent,interpolation='bilinear',origin='lower',cmap=CMAP,aspect=1.0)
    im.set_clim(-vmax,vmax)

    
    cb = plt.colorbar(orientation='vertical', shrink=0.785, pad = 0.007) 
    cb.set_ticks( np.linspace(-vmax,vmax,5) )
    cb.set_ticklabels( np.linspace(-vmax,vmax,5) )
    cb.ax.tick_params(axis='y', which='major', labelsize=FS-2)
    

#---------------------------------------------------------------------------------------
#   impose axis limits, add tick marks 
#--------------------------------------------------------------------------------------- 
    axis_lims = [0.,Lx/1000., 0., Ly/1000.]  
    plt.axis(axis_lims)
    
    majorLocator = MultipleLocator(100)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(25)
    axes.xaxis.set_major_locator(majorLocator)
    axes.xaxis.set_major_formatter(majorFormatter)
    # for the minor ticks, use no labels; default NullFormatter
    axes.xaxis.set_minor_locator(minorLocator)
    
    majorLocator = MultipleLocator(100)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(25)
    axes.yaxis.set_major_locator(majorLocator)
    axes.yaxis.set_major_formatter(majorFormatter)
    # for the minor ticks, use no labels; default NullFormatter
    axes.yaxis.set_minor_locator(minorLocator)

    axes.set_xlabel(r'$x\, \,km$',fontsize=FS)
    axes.set_ylabel(r'$y\, \,km$',fontsize=FS)
    title_string = r"$\zeta$      depth=%.2f m             t=%.2f  inertial periods" %( (Lz-z),t/iper)  
    axes.set_title(title_string, fontsize=FS-2)


#---------------------------------------------------------------------------------------
#  possibly, print the figure as an eps file
#---------------------------------------------------------------------------------------
    if do_print_eps == 1: 
        fname = plot_dir + 'zeta_' + str(islice).zfill(4) + '.eps' 
        plt.savefig(fname,dpi=600,bb_inches='tight') # save an eps w/ nice resolution  
        
#---------------------------------------------------------------------------------------
#  possibly, print the figure as a png file
#---------------------------------------------------------------------------------------
    if do_print_png == 1: 
        fname = plot_dir + 'zeta_' + str(islice).zfill(4) + '.png' 
        plt.savefig(fname,dpi=600,bb_inches='tight') # save a png  
        
    
    plt.close(fig)
    del axes
    plt.close('all')
    sys.exit()



