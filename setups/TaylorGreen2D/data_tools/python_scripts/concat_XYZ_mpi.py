#!
"""
Routine to concatenate XYZ netcdf files (3d)
in parallel using mpi4py.
usage example on morgon:  /usr/bin/mpirun -np 2 python concat_XYZ_mpi.py < /dev/null &> concat_log &
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
fnrs = sys.argv[7]               # filename root string

fnrs = fnrs + "_"                # more convenient to have the underscore here

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
	print(" concat_XYZ_mpi.py running on ",numprocs," processors on ", name) 
	print('!--------------------------------------------------------------------------!')


#-----------------------------------------------------------------------------------------
#  data/slice  location parameters
#-----------------------------------------------------------------------------------------
data_dir = root_dir + 'output/3D/'
slice_dir = root_dir + 'output/slices/3D/'
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
len = np.float( len(slices) )                             # length of global slice list, as float
local_len = np.int(np.ceil(len/numprocs))                 # nominal length of slice list on each processor
i0 = myid*local_len                                       # starting index for this processors slices
i1 = i0 + local_len                                       # (noninclusive) ending index
i1 = min(i1,last_slice+1)                                 # keep last list(s) from overshooting the end
my_slices = slices[i0:i1]                                 # list of slices assigned to processor myid

#print "myid = ",myid,"  ",my_slices                      # keep a record of the slice distribution...
comm.Barrier()                                            # all start together...


for islice in my_slices:
    print("......  myid, time slice: ",myid,"  ",islice)

#---------------------------------------------------------------------------------------
#  extract slice and concatenate files to create an XYZ snapshot  
#---------------------------------------------------------------------------------------
    # concatenate all the distributed 3d files...    
    cmd = "perl " + perl_dir +  "create_global_snapshot.pl 0 " + str(p1) + "  " + str(p2) + " " + data_dir[0:-1] + " XYZ_" + str(islice).zfill(6) + " " + slice_dir[0:-1] +  " " + fnrs + str(islice).zfill(6) + " " + str(myid)
    #print(cmd)
    os.system(cmd) 
    # condense the output filename and rename the file for convenience   
    cmd = "mv " + slice_dir + "XYZ_" + str(islice).zfill(6) +"_000000.nc " + slice_dir + fnrs + str(islice).zfill(6) + ".nc"
    os.system(cmd)
    
    # adjust storage order so that global file has same structure as flow_solve distributed files
    fname = slice_dir + fnrs + str(islice).zfill(6) + ".nc"
    cmd = "ncpdq -O -a timedimension,kdimension,jdimension,idimension " + fname + " " + fname
    os.system(cmd)
    
# all processors sync up before exiting routine    
comm.Barrier()
