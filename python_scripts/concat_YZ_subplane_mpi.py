#!
"""
Routine to concatenate EAST/WEST YZ netcdf files (2d) in parallel using mpi4py.
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
	print ("concat_YZ_subplane_mpi.py running on  ", numprocs, " processors on ", name) 
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
comm.Barrier()                                            # all start together...


#-----------------------------------------------------------------------------------------
#   flow_solve decomposition
#-----------------------------------------------------------------------------------------
locnx = int((nx-1)/p1)                      # x grid points per iprocessor, except p1 processor has 1 more
locnz = int((nz-1)/p2)                      # z grid points per jprocessor, except p2 processor has 1 more

jproc_start = int( np.floor( ( kstar[0] )/locnz ) )
jproc_end = int( np.floor( ( kstar[1] )/locnz ) )

iproc_east = int( np.floor( (istar[0])/locnx ) )      # the i processor of the east YZ planes
iproc_west = int( np.floor( (istar[1])/locnx ) ) # the i processor of the west YZ planes

if myid==0:
	print( '...         locnx,locnz             ',locnx,locnz )
	print( '...         jproc_start, jproc_end  ',jproc_start,jproc_end )
	print( '...         iproc_east, iproc_west  ',iproc_east, iproc_west )

for islice in my_slices:
    if( np.mod(islice,16)==0 ): print("......  myid, time slice: ",myid,"  ",islice)

#---------------------------------------------------------------------------------------
#  extract slice and concatenate files to create a YZplane snapshot 
#---------------------------------------------------------------------------------------
    cmd = "perl " + perl_dir + "concat_YZ_subplane.pl " + str(islice) + " " + str(iproc_east) + " "  + str(p2) + " " + data_dir[0:-1] + " east " + slice_dir[0:-1] +  " east " + str(myid) + " " + str(jproc_start) + " " + str(jproc_end)
    #print(cmd)
    os.system(cmd)
    
    cmd = "perl " + perl_dir + "concat_YZ_subplane.pl " + str(islice) + " " + str(iproc_west) + " "  + str(p2) + " " + data_dir[0:-1] + " west " + slice_dir[0:-1] +  " west " + str(myid) + " " + str(jproc_start) + " " + str(jproc_end)
    #print(cmd)
    os.system(cmd)

comm.Barrier() 
if myid==0:
	#---------------------------------------------------------------------------------------
	#  concatenate across time 
	#---------------------------------------------------------------------------------------
	cmd = "ncrcat -O " + slice_dir + "east_*.nc " + slice_dir + "east.nc"
	os.system(cmd)
	cmd = "ncrcat -O " + slice_dir + "west_*.nc " + slice_dir + "west.nc"
	os.system(cmd)

	#---------------------------------------------------------------------------------------
	#  clean up 
	#---------------------------------------------------------------------------------------
	cmd = "rm -f " + slice_dir + "east_*.nc"
	os.system(cmd)
	cmd = "rm -f " + slice_dir + "west_*.nc"
	os.system(cmd)
	
	cmd = "ls -lh " + slice_dir + "east.nc"
	os.system(cmd)
	cmd = "ls -lh " + slice_dir + "west.nc"
	os.system(cmd)

