#!
"""
Routine to concatenate XY netcdf files (2d)
in parallel using mpi4py.
usage example on morgon:  /usr/bin/mpirun -np 2 python concat_XY_mpi.py < /dev/null &> concat_log &
"""
#---------------------------------------------------------------------------------------
#  import and name the various modules I'll use
#---------------------------------------------------------------------------------------
import os,math,sys              # Import the os and math modules
import numpy as np
from scipy.io import netcdf
from mpi4py import MPI
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------
#  parse command line arguments
#-----------------------------------------------------------------------
script_name = sys.argv[0]
root_dir = sys.argv[1]           #  output needs to beneath root_dir
p1 = int(sys.argv[2])            # flow_solve decomposition parameter
p2 = int(sys.argv[3])            # flow_solve decomposition parameter
start_slice = int(sys.argv[4])
end_slice = int(sys.argv[5])
inc = int(sys.argv[6])
jproc = int(sys.argv[7])         # the jproc of the XY plane
fnrs = sys.argv[8]               # filename root string

slices = np.arange(start_slice,end_slice,inc)    #  list of time slices to process 


#-----------------------------------------------------------------------------------------
#   MPI initialization
#-----------------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
numprocs = MPI.COMM_WORLD.Get_size()    # number of processors executing this script
myid = MPI.COMM_WORLD.Get_rank()        # the processor id of this processor
name = MPI.Get_processor_name()         # myid's actual processor name, i.e. amigo.ucsd.edu
if myid==0:
	print('!--------------------------------------------------------------------------!') 
	print (" concat_XY_mpi.py running on  ",numprocs," processors on " ,name) 
	print('!--------------------------------------------------------------------------!')
  
#-----------------------------------------------------------------------------------------
#  data/slice  location parameters
#-----------------------------------------------------------------------------------------
data_dir = root_dir + 'output/2D/'
slice_dir = root_dir + 'output/slices/2D/'
perl_dir = root_dir + 'input/data_tools/perl_scripts/'
if myid==0:
    if not os.path.exists(slice_dir):
        cmd = 'mkdir -p ' + slice_dir
        os.system(cmd)
comm.Barrier()


#-----------------------------------------------------------------------------------------
#   divide up the slices across the numprocs processors
#   each mpi task makes its own list of slices: my_slices
#-----------------------------------------------------------------------------------------
last_slice=slices[-1]
len = float( len(slices) )                             # length of global slice list, as float
local_len = int(np.ceil(len/numprocs))                 # nominal length of slice list on each processor
i0 = myid*local_len                                       # starting index for this processors slices
i1 = i0 + local_len                                       # (noninclusive) ending index
i1 = min(i1,last_slice+1)                                 # keep last list(s) from overshooting the end
my_slices = slices[i0:i1]                                 # list of slices assigned to processor myid

#print "myid = ",myid,"  ",my_slices                      # keep a record of the slice distribution...
comm.Barrier()                                            # all start together...


for islice in my_slices:
    print("......  myid, time slice: ",myid,"  ",islice)

#---------------------------------------------------------------------------------------
#  extract slice and concatenate files to create a XYplane snapshot  p2=2 iproc=0
# ./concat_XY.pl tslice p1 jproc data_dir data_root out_dir out_root myid
#  ./concat_XY.pl 1050 2 1 '../../../output/2D' 'XYplane' '../../../output/slices/2D' 'XY' 0
#---------------------------------------------------------------------------------------
    cmd = "perl " + perl_dir +  "concat_XY.pl " + str(islice) + " " + str(p1) + " "  + str(jproc) + " " + data_dir[0:-1] + " " + fnrs + " " + slice_dir[0:-1] +  " " + fnrs + " " + str(myid)
    #print(cmd)
    os.system(cmd)    
    
# all processors sync up before exiting routine    
comm.Barrier()    
