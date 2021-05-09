#!
"""
Routine to concatenate YZ netcdf files (2d)
in parallel using mpi4py.
usage examples 
amigo :  mpirun -np 2 python concat_YZ_mpi.py
morgon:  /usr/bin/mpirun -np 2 python concat_YZ_mpi.py < /dev/null &> concat_log &
comet:  via run_concat.sh
"""
#---------------------------------------------------------------------------------------
#  import and name the various modules I'll use
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
from scipy.io import netcdf
from mpi4py import MPI
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#   MPI initialization
#-----------------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
numprocs = MPI.COMM_WORLD.Get_size()    # number of processors executing this script
myid = MPI.COMM_WORLD.Get_rank()        # the processor id of this processor
name = MPI.Get_processor_name()         # myid's actual processor name, i.e. amigo.ucsd.edu
if myid==0: 
 print " concat_YZ_mpi.py running on  %d processors on %s.\n" % (numprocs,name) 
  
#-----------------------------------------------------------------------------------------
#  data/slice  location parameters
#-----------------------------------------------------------------------------------------
slice_dir = '../../output/slices/2D/'
data_dir = '../../output/2D/'
if myid==0:
    if not os.path.exists(slice_dir):
        cmd = 'mkdir -p ' + slice_dir
        os.system(cmd)
comm.Barrier()

slices = np.arange(0,176,1)                   #  list of time slices to process 
last_slice=slices[-1]

#-----------------------------------------------------------------------------------------
#   divide up the slices across the numprocs processors
#   each mpi task makes its own list of slices: my_slices
#-----------------------------------------------------------------------------------------
len = np.float( len(slices) )                             # length of global slice list, as float
local_len = np.int(np.ceil(len/numprocs))                 # nominal length of slice list on each processor
i0 = myid*local_len                                       # starting index for this processors slices
i1 = i0 + local_len                                       # (noninclusive) ending index
i1 = min(i1,last_slice+1)                                 # keep last list(s) from overshooting the end
my_slices = slices[i0:i1]                                 # list of slices assigned to processor myid

#print "myid = ",myid,"  ",my_slices                       # keep a record of the slice distribution...
comm.Barrier()                                            # all start together...


#-----------------------------------------------------------------------------------------
#   flow_solve decomposition
#-----------------------------------------------------------------------------------------
p1 = 8        #  number of processors splitting nx and ny  
p2 = 24       #  number of processors splitting ny and nz   
nx=256 ; ny=313; nz=265                #  number of grid points
istar=[102,154]                        # i limits specified in io_params
jstar=[129,182]                        # j limits specified in io_params
kstar=[1,265]                          # j limits specified in io_params
locnx = nx/p1                          # x grid points per iprocessor
locnz = (nz-1)/p2                      # z grid points per jprocessor, except p2 processor has 1 more

jproc_start = np.int( np.floor( ( kstar[0]-1.e-6 )/locnz ) )
jproc_end = np.int( np.floor( ( kstar[1]-1.e-6 )/locnz ) )

for islice in my_slices:
    print "processing time slice: %d %d" % (myid,islice)

#---------------------------------------------------------------------------------------
#  extract slice and concatenate files to create a YZplane snapshot 
#---------------------------------------------------------------------------------------
    iproc = np.int( np.floor( (istar[0]-1.e-6)/locnx ) )  # the i processor of the X=X0 YZ planes
    cmd = "perl ./concat_YZ_subplane.pl " + str(islice) + " " + str(iproc) + " "  + str(p2) + " " + data_dir[0:-1] + " YZ_X=X0 " + slice_dir[0:-1] +  " YZ_X=X0 " + str(myid) + " " + str(jproc_start) + " " + str(jproc_end)
    print cmd
    os.system(cmd)
    
    iproc = np.int( np.floor( (istar[1]-1.e-6)/locnx ) )  # the i processor of the X=X1 YZ planes
    cmd = "perl ./concat_YZ_subplane.pl " + str(islice) + " " + str(iproc) + " "  + str(p2) + " " + data_dir[0:-1] + " YZ_X=X1 " + slice_dir[0:-1] +  " YZ_X=X1 " + str(myid) + " " + str(jproc_start) + " " + str(jproc_end)
    print cmd
    os.system(cmd)
