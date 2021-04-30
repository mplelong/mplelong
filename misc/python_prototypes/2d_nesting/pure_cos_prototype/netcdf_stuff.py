def write_netcdf_2d(arrays,xvec,zvec,tval,filename,varnames,dims):
#
#	Simple routine to write 2d spatial data to the netcdf file filename: 
#
#     arrays is a list of 2d data arrays
#     2d data stored in python arrays with layout f[x,z] with names varnames and dimensions dims
#     The 1d arrays x,z are assumed to be spatial arrays in [m]
#
#
#   5 flavors of netcdf files:
#     (NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET, NETCDF3_64BIT_DATA, NETCDF4_CLASSIC, and NETCDF4)
#     Default NETCDF4 files use the version 4 disk format (HDF5) and use the new features of the version 4 API.
#     To see how a given file is formatted, you can examine the data_model attribute.
#
#   Modes of opening netcdf files
#     r means read-only; no data can be modified. 
#     w means write; a new file is created, an existing file with the same name is deleted. 
#     a and r+ mean append (in analogy with serial files); an existing file is opened for reading and writing
#
# 
#------------------------------------------------------------------------------------------------------------------
	import netCDF4 as nc
	import numpy as np

	success = False    # switch to True if all steps are completed
	
	fmt = "NETCDF4"
	mode = 'w'
	create_dims = True
	dtype = np.float64
	nvars = len(arrays)
			
	nx = np.size(xvec)
	nz = np.size(zvec)
	nt = 1
	
	# open or create the netcdf file with the appropriate name and mode	
	nc = nc.Dataset(filename, mode, format=fmt)   

	if( create_dims ):
		# Every dimension has a name and length. If we set the dimension length to be 0 or None, 
		# then it takes it as of unlimited size and can grow.
		tdim = nc.createDimension("tdim", nt)    # None triggers unlimited dimension so I can concatenate files across time later
		idim = nc.createDimension("idim", nx)
		kdim = nc.createDimension("kdim", nz)
	
		# A variable has a name, a type, a shape and some data values. The shape of the variable can be stated 	
		# using the tuple of the dimension names. Can also add attributes.
		t = nc.createVariable('t', dtype, ("tdim",))
		x = nc.createVariable('x', dtype, ("idim",))
		z = nc.createVariable('z', dtype, ("kdim",))
		t.units = "s"
		x.units = "m"
		z.units = "m"
		
		# Write the dimension data
		t[:] = tval
		x[:] = xvec[:]
		z[:] = zvec[:]
		
	
	# create the data variable, give it a units attribute and write the data
	for i in range(nvars):
		data = nc.createVariable(varnames[i], dtype, ("tdim","kdim","idim",))
		data.units = dims[i]
		data[:,:] = np.transpose(arrays[i])
	
	
	# print the Dataset object to see what we've got and then close the file
	# print(nc)
	nc.close()
	
	print(filename,' has been created, written to and closed.')
	success = True
	
	return success

def read_netcdf_2d_coords(filename,c1,c2):
#
#	Read the coordinates x,z and t from filename
#
	import netCDF4 as nc
	import numpy as np	
	ds = nc.Dataset(filename)
	#t = ds['time'][:]
	t = ds['t'][:]
	x = ds[c1][:]
	z = ds[c2][:]
	return t,x,z   

def read_netcdf_2d_var(varname,filename,tslice):
#
#	Read a 2d array f[t,z,x] from an netcdf file, return the array f[z,x].
#
	import netCDF4 as nc
	import numpy as np	
	ds = nc.Dataset(filename)
	f = ds[varname][tslice,:,:]
	return np.squeeze(f) 
	
def read_netcdf_2d_x0(varname,filename,tslice,xindex):
#
#	Read a z profile f[t,z,x0] from an netcdf file, return the array f[z].
#
	import netCDF4 as nc
	import numpy as np	
	ds = nc.Dataset(filename)
	f = ds[varname][tslice,:,xindex]
	return np.squeeze(f)  
	
def read_netcdf_2d_z0(varname,filename,tslice,zindex):
#
#	Read a z profile f[t,z0,:] from an netcdf file, return the array f[x].
#
	import netCDF4 as nc
	import numpy as np	
	ds = nc.Dataset(filename)
	f = ds[varname][tslice,zindex,:]
	return np.squeeze(f) 

