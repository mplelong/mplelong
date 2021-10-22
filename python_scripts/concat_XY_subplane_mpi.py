#!
"""
Routine to concatenate BOTTOM/TOP XY netcdf files (2d) in parallel using mpi4py.
"""
#---------------------------------------------------------------------------------------
#  import and name the various modules I'll use
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
from scipy.io import netcdf
from mpi4py import MPI
from data_processing_utilities import parse_problem_params
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------
#  parse command line arguments
#-----------------------------------------------------------------------
script_name = sys.argv[0]
root_dir = sys.argv[1]           #  output needs to beneath root_dir
root_dir = root_dir + '/'
start_slice = int(sys.argv[2])
end_slice = int(sys.argv[3])
inc = int(sys.argv[4])
i0 = int(sys.argv[5])
i1 = int(sys.argv[6])
j0 = int(sys.argv[7])
j1 = int(sys.argv[8])
k0 = int(sys.argv[9])
k1 = int(sys.argv[10])

slice_dir = root_dir + 'output/slices/2D/'
data_dir = root_dir + 'output/2D/'
perl_dir = root_dir + 'input/data_tools/perl_scripts/'

slices = np.arange(start_slice,end_slice,inc)    #  list of time slices to process 
last_slice=slices[-1]

istar=[i0,i1]   # i limits specified in io_params, read from command line
jstar=[j0,j1]   # j limits specified in io_params
kstar=[k0,k1]   # k limits specified in io_params

problem_params = parse_problem_params(root_dir)
[runlabel,restart_flag,do_second_scalar,AB_order,p1,p2,nx,ny,nz,dt,t0,tf,Lx,Ly,Lz, \
x_periodic,y_periodic,z_periodic,s1_def,s2_def,user_forcing_flag,rho0,g,f0,nu,     \
kappa1,kappa2,high_order_flag,p,T_diff] = problem_params


#-----------------------------------------------------------------------------------------
#   MPI initialization
#-----------------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
numprocs = MPI.COMM_WORLD.Get_size()    # number of processors executing this script
myid = MPI.COMM_WORLD.Get_rank()        # the processor id of this processor
name = MPI.Get_processor_name()         # myid's actual processor name, i.e. amigo.ucsd.edu
if myid==0:
	print('!--------------------------------------------------------------------------!') 
	print ("concat_XY_subplane_mpi.py running on  ", numprocs, " processors on ", name) 
	print('!--------------------------------------------------------------------------!')
	print('...         flow_solve root directory:', root_dir)
	print('...         input data directory     :', data_dir)
	print('...         output slice directory   :', slice_dir)
	print('...         flow_solve decomposition :', p1,p2)
	print('...         flow_solve resolution    :', nx,ny,nz)
	print('...         time slice parameters    :', start_slice,end_slice,inc)
	print('...         child grid indices       :', istar,jstar,kstar) 
  
if myid==0:
    if not os.path.exists(slice_dir):
        cmd = 'mkdir -p ' + slice_dir
        os.system(cmd)
comm.Barrier()

#-----------------------------------------------------------------------------------------
#   divide up the slices across the numprocs processors
#   each mpi task makes its own list of slices: my_slices
#-----------------------------------------------------------------------------------------
len = float( len(slices) )                                # length of global slice list, as float
local_len = int(np.ceil(len/numprocs))                    # nominal length of slice list on each processor
i0 = myid*local_len                                       # starting index for this processors slices
i1 = i0 + local_len                                       # (noninclusive) ending index
i1 = min(i1,last_slice+1)                                 # keep last list(s) from overshooting the end
my_slices = slices[i0:i1]                                 # list of slices assigned to processor myid

#print "myid = ",myid,"  ",my_slices                       # keep a record of the slice distribution...
comm.Barrier()


#-----------------------------------------------------------------------------------------
#   flow_solve decomposition etc
#-----------------------------------------------------------------------------------------
locnx = int((nx-1)/p1)                      # x grid points per iprocessor, except p1 processor has 1 more
locnz = int((nz-1)/p2)                      # z grid points per jprocessor, except p2 processor has 1 more


iproc_start = int( np.floor( ( istar[0] )/locnx ) )
iproc_end = int( np.floor( ( istar[1] )/locnx ) )

jproc_bot = int( np.floor( (kstar[0])/locnz ) )                   # the j processor of the bottom XY plane
jproc_top = min( [int( np.floor( (kstar[1])/locnz ) ),p2-1] )  # the j processor of the top XY plane

if myid==0:
	print( '...         locnx,locnz             ',locnx,locnz )
	print( '...         jproc_bot, jproc_top    ',jproc_bot,jproc_top )
	print( '...         iproc_start, iproc_end  ',iproc_start, iproc_end )
	

for islice in my_slices:
    if( np.mod(islice,16)==0 ): print("......  myid, time slice: ",myid,"  ",islice)

#---------------------------------------------------------------------------------------
#  extract slice and concatenate files to create a XYplane snapshot  p2=2 iproc=0
# ./concat_XY.pl tslice p1 jproc data_dir data_root out_dir out_root myid iproc_start iproc_end
#---------------------------------------------------------------------------------------
    cmd = "perl " + perl_dir + "concat_XY_subplane.pl " + str(islice) + " " + str(p1) + " "  + str(jproc_bot) + " " + data_dir[0:-1] + " bottom " + slice_dir[0:-1] +  " bottom " + str(myid) + " " + str(iproc_start) + " " + str(iproc_end)
    #print(cmd)
    os.system(cmd)
    
    cmd = "perl " + perl_dir + "concat_XY_subplane.pl " + str(islice) + " " + str(p1) + " "  + str(jproc_top) + " " + data_dir[0:-1] + " top " + slice_dir[0:-1] +  " top " + str(myid) + " " + str(iproc_start) + " " + str(iproc_end)
    #print(cmd)
    os.system(cmd)


comm.Barrier() 
if myid==0:
	#---------------------------------------------------------------------------------------
	#  concatenate across time 
	#---------------------------------------------------------------------------------------
	cmd = "ncrcat -O " + slice_dir + "bottom_*.nc " + slice_dir + "bottom.nc"
	os.system(cmd)
	cmd = "ncrcat -O " + slice_dir + "top_*.nc " + slice_dir + "top.nc"
	os.system(cmd)

	#---------------------------------------------------------------------------------------
	#  clean up 
	#---------------------------------------------------------------------------------------
	cmd = "rm -f " + slice_dir + "bottom_*.nc"
	os.system(cmd)
	cmd = "rm -f " + slice_dir + "top_*.nc"
	os.system(cmd)
	
	cmd = "ls -lh " + slice_dir + "bottom.nc"
	os.system(cmd)
	cmd = "ls -lh " + slice_dir + "top.nc"
	os.system(cmd)    
