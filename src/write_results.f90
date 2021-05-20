subroutine write_results
	use io_params 
	use etc,           only: istep,istart 
	use mpi_params,    only: myid,comm,ierr
  
	implicit none
	integer               :: fid
	integer,save          :: jstart(maxsets)
	logical,save          :: first_entry=.TRUE.
 
	if( first_entry ) then
		jstart(:)=istart      ! generally, all output sets start up at istep=istart
		first_entry=.false.
	endif
   
	!-------------------------------------------------------
	! loop through all the user-specified file sets...
	!-------------------------------------------------------
	do fid=1,num_file_sets 
		if( istep >= jstart(fid) ) then
			if( mod( (istep-jstart(fid)),nsteps(fid) ) == 0    &
			   .or. istep == jstart(fid) ) then
    
				!-------------------------------------------------------
    			! initialize netcdf file if necessary
    			!-------------------------------------------------------
    			if( trim(mode(fid)) == 'new' .or. istep == jstart(fid)) then      
     				call init_netcdf(fid) 
    			endif 
   
    			!-------------------------------------------------------
    			! file already initialized, just write current data
    			!-------------------------------------------------------
    			call write_netcdf(fid)
   
			endif   
		endif  
	enddo   ! end loop through file sets

	call mpi_barrier(comm,ierr)
 return
end subroutine write_results



  
  
subroutine init_netcdf(fid)
	use io_params
	use decomposition_params
	use methods_params,         only: do_second_scalar
	use dependent_variables
	use independent_variables
	use etc,                    only: istep
	use mpi_params,             only: myid,comm,ierr
      
	implicit none
	include 'netcdf.inc'

	character(len=80)              :: filename
	character(len=3)               :: dimstr
	character(len=80),save         :: topdir='output/'
	character(len=80)              :: s1_name,s2_name
	character(len=80)              :: s1_units,s2_units
	integer                        :: i,j,k
	integer                        :: fid,ncid,rcode
	integer                        :: xVarID,iid
	integer                        :: yVarID,jid
	integer                        :: zVarID,kid
	integer                        :: tVarID,timeid
	integer                        :: All_nD_VarSIZE(4),ndims
	integer                        :: uVarID,vVarID,wVarID 
	integer                        :: ustarVarID,vstarVarID,wstarVarID
	integer                        :: s1VarID,s2VarID
	integer                        :: s1_barVarID,s2_barVarID
	integer                        :: divustarVarID,phiVarID,pdVarID
	integer                        :: uderivVarID,vderivVarID,wderivVarID
	integer                        :: s1derivVarID,s2derivVarID
	integer                        :: s1_name_len,s2_name_len
	integer                        :: s1_units_len,s2_units_len
	integer                        :: npts,start1D(1),count1D(1)
	integer                        :: counter,offset
	real(kind=8),allocatable,save  :: scratch(:)
	logical,save                   :: first_entry=.TRUE.

	if( first_entry ) then
		npts = maxval((/nx,ny,nz/))
		allocate( scratch(npts)  )
		scratch = 0.d0
		first_entry=.FALSE.
	endif
 
	call process_dimensions(fid,dimstr)  !! figure out if file should be 1D,2D,3D etc
 
	call construct_filename(fid,dimstr,topdir,istep,fullname(fid))
  
	!---------------------------------------------------------------
	! specify the start and count arrays used in netcdf write calls
	! set the logical variable do_write(fid), this routine will assume
	! that there is data to write for myid and detect if this is
	! really the case. If not, it will set do_write(fid) to .FALSE.
	!---------------------------------------------------------------
	call set_local_indices_count(fid)
	time_counter(fid)=1    ! always 1 for initialization
	if( .NOT. do_write(fid) ) goto 999
 
 
	!--------------------------------------------------------------- 
	!  Open a netcdf file
	!---------------------------------------------------------------
	rcode=NF_CREATE(trim(fullname(fid)),NF_NOCLOBBER,ncid)
	if(rcode.ne.NF_NOERR) then
		write(0,*) '... ERROR OPENING NETCDF FILE: init_netcdf ',trim(fullname(fid))
		write(0,*) '... myid, rcode ',myid, nf_strerror(rcode)
		stop
	endif
  

	!--------------------------------------------------------------- 
	! Define (and order) the spatial dimensions
	!--------------------------------------------------------------- 
	rcode=NF_DEF_DIM(ncid,'idimension',count(dimid_x(fid),fid),iid)
	if (rcode.ne.NF_NOERR) write(0,*) myid,  &
         ': NetCDF Error: NF_DEF_DIM: idimension', rcode
	All_nD_VarSIZE(dimid_x(fid))=iid

	rcode=NF_DEF_DIM(ncid,'jdimension',count(dimid_y(fid),fid),jid)
	if (rcode.ne.NF_NOERR) write(0,*) myid,  &
          ': NetCDF Error: NF_DEF_DIM: jdimension', rcode
	All_nD_VarSIZE(dimid_y(fid))=jid

	rcode=NF_DEF_DIM(ncid,'kdimension',count(dimid_z(fid),fid),kid)
	if (rcode.ne.NF_NOERR) write(0,*) myid,   &
          ': NetCDF Error: NF_DEF_DIM: kdimension', rcode
	All_nD_VarSIZE(dimid_z(fid))=kid
   
	rcode=NF_DEF_DIM(ncid,'timedimension',NF_UNLIMITED,timeid)
	if (rcode.ne.NF_NOERR) write(0,*) myid,   &
          ': NetCDF Error: NF_DEF_DIM: timedimension', rcode
	All_nD_VarSIZE(dimid_t(fid))=timeid
	ndims = nspace(fid) + 1
  

	!---------------------------------------------------------------
	!  Define x,y,z grid position  and time variables
	!---------------------------------------------------------------

	!!**X*****
	rcode=NF_DEF_VAR(ncid,'x',NF_DOUBLE,1,iid,xVarID)
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: x', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,xVarID,'long_name',1,'x')
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: x', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,xVarID,'units',1,'m')
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: x', rcode

	!!**Y*****
	rcode=NF_DEF_VAR(ncid,'y',NF_DOUBLE,1,jid,yVarID)
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: y', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,yVarID,'long_name',1,'y')
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: y', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,yVarID,'units',1,'m')
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: y', rcode
           
	!!**Z*****
	rcode=NF_DEF_VAR(ncid,'z',NF_DOUBLE,1,kid,zVarID)
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: z', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,zVarID,'long_name',1,'z')
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: z', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,zVarID,'units',1,'m')
	if (rcode.ne.NF_NOERR) print *,myid,   &
       ': NetCDF Error: NF_PUT_ATT_TEXT: z', rcode

	!!**time*****
	rcode=NF_DEF_VAR(ncid,'time',NF_DOUBLE,1,timeid,tVarID)
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_DEF_VAR: time', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,tVarID,'long_name',4,'time')
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: time', rcode
	rcode=NF_PUT_ATT_TEXT(ncid,tVarID,'units',7,'seconds')
	if (rcode.ne.NF_NOERR) print *,myid,   &
      ': NetCDF Error: NF_PUT_ATT_TEXT: time', rcode     

  
	!-----------------------------------------------
	! Set the appropriate labels for s1 and s2.
	!-----------------------------------------------
	if(trim(scalar_kind(1)) == 't') then
		s1_name = 'Temperature'
		s1_name_len = 11
  	  s1_units = 'deg C'
  	  s1_units_len = 5
 	 endif
 	 if(trim(scalar_kind(2)) == 't') then
 	   s2_name = 'Temperature'
  	  s2_name_len = 11
  	  s2_units = 'deg C'
  	  s2_units_len = 5
 	 endif

  	if(trim(scalar_kind(1)) == 's') then
  	  s1_name = 'Salinity'
  	  s1_name_len = 8
  	  s1_units = 'psu'
  	  s1_units_len = 3
 	 endif
 	 if(trim(scalar_kind(2)) == 's') then
 	   s2_name = 'Salinity'
 	   s2_name_len = 8
 	   s2_units = 'psu'
 	   s2_units_len = 3
 	 endif
       
 	 if(trim(scalar_kind(1)) == 'p') then
  	  s1_name = 'Passive Tracer'
  	  s1_name_len = 14
  	  s1_units = 'Concentration'
  	  s1_units_len = 13
 	 endif
 	 if(trim(scalar_kind(2)) == 'p') then
  	  s2_name = 'Passive Tracer'
  	  s2_name_len = 14
  	  s2_units = 'Concentration'
  	  s2_units_len = 13
 	 endif
             
 	 if(trim(scalar_kind(1)) == 'r') then
 	   s1_name = 'Density'
 	   s1_name_len = 7
 	   s1_units = 'kg/m3'
 	   s1_units_len = 5
 	 endif
 	 if(trim(scalar_kind(2)) == 'r') then
 	   s2_name = 'Density'
 	   s2_name_len = 7
 	   s2_units = 'kg/m3'
 	   s2_units_len = 5
	  endif
  
 	 if(trim(scalar_kind(2)) == 'x') then
 	   s2_name = 'Sediment Concentration'
 	   s2_name_len = 22
 	   s2_units = '[1]'
 	   s2_units_len = 3
	  endif
      

	!!
	! scalar (1)
	!!   
    if( variable_key(4,fid) /= 0 ) then    
    	rcode=NF_DEF_VAR(ncid,'s1_bar',NF_DOUBLE,1,kid,s1_barVarID)
      	if (rcode.ne.NF_NOERR) print *,myid,   &
          ': NetCDF Error: NF_DEF_VAR: s1_bar', rcode       
     	 rcode=NF_DEF_VAR(ncid,'s1',NF_DOUBLE,ndims,All_nD_VarSIZE(1:ndims),s1VarID)
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: s1', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,s1VarID,'long_name',s1_name_len,s1_name)
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s1', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,s1VarID,'units',s1_units_len,s1_units)
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s1', rcode
    endif

	!!
	!scalar (2)
	!!
    if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
     	 rcode=NF_DEF_VAR(ncid,'s2_bar',NF_DOUBLE,1,kid,s2_barVarID)
     	 if (rcode.ne.NF_NOERR) print *,myid,   &
          ': NetCDF Error: NF_DEF_VAR: s2_bar', rcode
     	 rcode=NF_DEF_VAR(ncid,'s2',NF_DOUBLE,ndims,All_nD_VarSIZE,s2VarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: s2', rcode           
     	 rcode=NF_PUT_ATT_TEXT(ncid,s2VarID,'long_name',s2_name_len,s2_name)                           
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s2', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,s2VarID,'units',s2_units_len,s2_units)
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: s2', rcode
    endif

	!!
	!x velocity = u
	!!
    if( variable_key(1,fid) == 1 ) then
      	rcode=NF_DEF_VAR(ncid,'u',NF_DOUBLE,ndims,All_nD_VarSIZE,uVarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: u', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,uVarID,'long_name',1,'u')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: u', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,uVarID,'units',3,'m/s')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: u', rcode
    endif
	!!
	!y velocity = v
	!! 
    if( variable_key(2,fid) == 1 ) then
    	  rcode=NF_DEF_VAR(ncid,'v',NF_DOUBLE,ndims,All_nD_VarSIZE,vVarID)
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: v', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,vVarID,'long_name',1,'v')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: v', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,vVarID,'units',3,'m/s')
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: v', rcode
    endif
	!!
	!z velocity = w
	!!
    if( variable_key(3,fid) == 1 ) then
      	rcode=NF_DEF_VAR(ncid,'w',NF_DOUBLE,ndims,All_nD_VarSIZE,wVarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: w', rcode
      	rcode=NF_PUT_ATT_TEXT(ncid,wVarID,'long_name',1,'w')
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: w', rcode
      	rcode=NF_PUT_ATT_TEXT(ncid,wVarID,'units',3,'m/s')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: w', rcode
    endif
        
	!!
	!divustar = divustar
	!! 
    if( variable_key(6,fid) == 1 ) then
     	 rcode=NF_DEF_VAR(ncid,'divustar',NF_DOUBLE,ndims,All_nD_VarSIZE,divustarVarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: divustar', rcode
      	rcode=NF_PUT_ATT_TEXT(ncid,divustarVarID,'long_name',9,'div ustar')                            
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: divustar', rcode
      	rcode=NF_PUT_ATT_TEXT(ncid,divustarVarID,'units',4,'1/s')
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: divustar', rcode
    endif
    
	!!
	!pressure variable = phi
	!! 
    if( variable_key(7,fid) == 1 ) then
      	rcode=NF_DEF_VAR(ncid,'phi',NF_DOUBLE,ndims,All_nD_VarSIZE,phiVarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: phi', rcode
      	rcode=NF_PUT_ATT_TEXT(ncid,phiVarID,'long_name',8,'pressure')                            
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: phi', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,phiVarID,'units',6,'kg/ms2')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: phi', rcode
    endif
    
    !!
	!x velocity = ustar
	!!
    if( variable_key(9,fid) == 1 ) then
      	rcode=NF_DEF_VAR(ncid,'ustar',NF_DOUBLE,ndims,All_nD_VarSIZE,ustarVarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: ustar', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,ustarVarID,'long_name',5,'ustar')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: ustar', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,ustarVarID,'units',3,'m/s')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: ustar', rcode
    endif
    
    !!
	!y velocity = vstar
	!!
    if( variable_key(10,fid) == 1 ) then
      	rcode=NF_DEF_VAR(ncid,'vstar',NF_DOUBLE,ndims,All_nD_VarSIZE,vstarVarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: vstar', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,vstarVarID,'long_name',5,'vstar')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: vstar', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,vstarVarID,'units',3,'m/s')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: vstar', rcode
    endif
    
    !!
	!z velocity = wstar
	!!
    if( variable_key(11,fid) == 1 ) then
      	rcode=NF_DEF_VAR(ncid,'wstar',NF_DOUBLE,ndims,All_nD_VarSIZE,wstarVarID)
      	 if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_DEF_VAR: wstar', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,wstarVarID,'long_name',5,'wstar')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: wstar', rcode
     	 rcode=NF_PUT_ATT_TEXT(ncid,wstarVarID,'units',3,'m/s')
     	  if (rcode.ne.NF_NOERR) print *,myid,   &
           ': NetCDF Error: NF_PUT_ATT_TEXT: wstar', rcode
    endif
    
    
    !----------------------------------------------------------
    ! child grid needs extra variables: 
    ! the normal derivatives at the 6 boundary faces
    !-----------------------------------------------------------
    if( t_secs >= t_start_child .and. do_child_grid .and. filename_root(fid) .NE. 'XYZ_child_' ) then
    
    	if( filename_root(fid)=='east' .or. filename_root(fid)=='west') then
    	
    		rcode=NF_DEF_VAR(ncid,'u_x',NF_DOUBLE,ndims,All_nD_VarSIZE,uderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,uderivVarID,'long_name',3,'u_x')
    		rcode=NF_PUT_ATT_TEXT(ncid,uderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'v_x',NF_DOUBLE,ndims,All_nD_VarSIZE,vderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,vderivVarID,'long_name',3,'v_x')
    		rcode=NF_PUT_ATT_TEXT(ncid,vderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'w_x',NF_DOUBLE,ndims,All_nD_VarSIZE,wderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,wderivVarID,'long_name',3,'w_x')
    		rcode=NF_PUT_ATT_TEXT(ncid,wderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'s1_x',NF_DOUBLE,ndims,All_nD_VarSIZE,s1derivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,s1derivVarID,'long_name',3,'s1_x')
    		rcode=NF_PUT_ATT_TEXT(ncid,s1derivVarID,'units',6,'[s1]/m]')
    		
    		if(write_s2_child_grid) then
    			rcode=NF_DEF_VAR(ncid,'s2_x',NF_DOUBLE,ndims,All_nD_VarSIZE,s2derivVarID)
    			rcode=NF_PUT_ATT_TEXT(ncid,s2derivVarID,'long_name',4,'s2_x')
    			rcode=NF_PUT_ATT_TEXT(ncid,s2derivVarID,'units',6,'[s2]/m]')
    		endif	
    			
    	elseif( filename_root(fid)=='south' .or. filename_root(fid)=='north') then
    
    		rcode=NF_DEF_VAR(ncid,'u_y',NF_DOUBLE,ndims,All_nD_VarSIZE,uderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,uderivVarID,'long_name',3,'u_y')
    		rcode=NF_PUT_ATT_TEXT(ncid,uderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'v_y',NF_DOUBLE,ndims,All_nD_VarSIZE,vderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,vderivVarID,'long_name',3,'v_y')
    		rcode=NF_PUT_ATT_TEXT(ncid,vderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'w_y',NF_DOUBLE,ndims,All_nD_VarSIZE,wderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,wderivVarID,'long_name',3,'w_y')
    		rcode=NF_PUT_ATT_TEXT(ncid,wderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'s1_y',NF_DOUBLE,ndims,All_nD_VarSIZE,s1derivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,s1derivVarID,'long_name',3,'s1_y')
    		rcode=NF_PUT_ATT_TEXT(ncid,s1derivVarID,'units',6,'[s1]/m]')
    		
    		if(write_s2_child_grid) then
    			rcode=NF_DEF_VAR(ncid,'s2_y',NF_DOUBLE,ndims,All_nD_VarSIZE,s2derivVarID)
    			rcode=NF_PUT_ATT_TEXT(ncid,s2derivVarID,'long_name',4,'s2_y')
    			rcode=NF_PUT_ATT_TEXT(ncid,s2derivVarID,'units',6,'[s2]/m]')
    		endif
    	
    	
    	elseif( filename_root(fid)=='bottom' .or. filename_root(fid)=='top') then
    	
    		rcode=NF_DEF_VAR(ncid,'u_z',NF_DOUBLE,ndims,All_nD_VarSIZE,uderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,uderivVarID,'long_name',3,'u_z')
    		rcode=NF_PUT_ATT_TEXT(ncid,uderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'v_z',NF_DOUBLE,ndims,All_nD_VarSIZE,vderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,vderivVarID,'long_name',3,'v_z')
    		rcode=NF_PUT_ATT_TEXT(ncid,vderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'w_z',NF_DOUBLE,ndims,All_nD_VarSIZE,wderivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,wderivVarID,'long_name',3,'w_z')
    		rcode=NF_PUT_ATT_TEXT(ncid,wderivVarID,'units',3,'1/s]')
    		
    		rcode=NF_DEF_VAR(ncid,'s1_z',NF_DOUBLE,ndims,All_nD_VarSIZE,s1derivVarID)
    		rcode=NF_PUT_ATT_TEXT(ncid,s1derivVarID,'long_name',3,'s1_z')
    		rcode=NF_PUT_ATT_TEXT(ncid,s1derivVarID,'units',6,'[s1]/m]')
    		
    		if(write_s2_child_grid) then
    			rcode=NF_DEF_VAR(ncid,'s2_z',NF_DOUBLE,ndims,All_nD_VarSIZE,s2derivVarID)
    			rcode=NF_PUT_ATT_TEXT(ncid,s2derivVarID,'long_name',4,'s2_z')
    			rcode=NF_PUT_ATT_TEXT(ncid,s2derivVarID,'units',6,'[s2]/m]')
    		endif
    		
    	endif
    	
    endif ! end do_child_grid block
    
    
	!------------------------------------------------------------------------
	!End define mode
	!------------------------------------------------------------------------
    rcode=NF_ENDDEF(ncid)
    if(rcode.ne.NF_NOERR) print *,myid,'ERROR  LEAVING DEFINE MODE'
	!------------------------------------------------------------------------
	!End define mode
	!------------------------------------------------------------------------

	!    store x values in temp array
    counter = 1
    offset = global_x_indices(START,YBLOCK,myid) - 1
    do i = my_x0(fid)+offset,my_x1(fid)+offset,my_xinc(fid)
    	scratch(counter) = x(i) ! [m]
    	counter = counter + 1
    enddo

	!    Put dimension data in netcdf file
    start1D(1) = 1
    count1D(1) = my_nx(fid)
    rcode=NF_PUT_VARA_DOUBLE(ncid,xVarID,start1D,count1D,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR xVarID: ',rcode


	!    store y values in temp array
    counter = 1
    offset = global_y_indices(START,YBLOCK,myid) - 1
    do j = my_y0(fid)+offset,my_y1(fid)+offset,my_yinc(fid)
    	scratch(counter) = y(j) ! [m]
    	counter = counter + 1
    enddo

	!    Put dimension data in netcdf file
    start1D(1) = 1
    count1D(1) = my_ny(fid)
    rcode=NF_PUT_VARA_DOUBLE(ncid,yVarID,start1D,count1D,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR yVarID: ',rcode,filename,count1D(1)

    

	!    store z values in temp array
    counter = 1
    offset = global_z_indices(START,YBLOCK,myid) - 1
    do k = my_z0(fid)+offset,my_z1(fid)+offset,my_zinc(fid)
    	scratch(counter) = z(k) ! [m]
    	counter = counter + 1
    enddo

	!    Put dimension data in netcdf file
    start1D(1) = 1
    count1D(1) = my_nz(fid)
    rcode=NF_PUT_VARA_DOUBLE(ncid,zVarID,start1D,count1D,scratch)
    if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR zVarID: ',rcode
      

	!     Write s1_bar if s1 itself is to be written     
    if( variable_key(4,fid) /= 0 ) then
    	counter = 1
    	offset = global_z_indices(START,YBLOCK,myid) - 1
    	do k = my_z0(fid)+offset,my_z1(fid)+offset,my_zinc(fid)
    		scratch(counter) = s1_bar(k,1)
    		counter = counter + 1
    	enddo
    	rcode=NF_PUT_VARA_DOUBLE(ncid,s1_barVarID,start1D,count1D,scratch)
    	if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR s1_barVarID: ',rcode
    endif
      
	!     Write s2_bar if s2 itself is to be written     
    if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
    	counter = 1
    	offset = global_z_indices(START,YBLOCK,myid) - 1
    	do k = my_z0(fid)+offset,my_z1(fid)+offset,my_zinc(fid)
    		scratch(counter) = s2_bar(k,1)
    		counter = counter + 1
    	enddo
    	rcode=NF_PUT_VARA_DOUBLE(ncid,s2_barVarID,start1D,count1D,scratch)
    	if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR s2_barVarID: ',rcode
    endif      
      
    rcode=NF_CLOSE(ncid)
    if(rcode.ne.NF_NOERR) print *,myid,'ERROR CLOSING NETCDF FILE'
        
999 continue
    call mpi_barrier(comm,ierr)  
 return
end subroutine init_netcdf




subroutine write_netcdf(fid)
	use decomposition_params
 	use dimensional_scales
 	use dependent_variables
 	use independent_variables
 	use intermediate_variables, only: div_u,phi,ustar,vstar,wstar
 	use methods_params,         only: do_second_scalar
 	use mpi_params
 	use io_params
 	use etc
      
 	implicit none
 	include 'netcdf.inc'
 

 	character(len=80)         :: err_msg
 	integer                   :: i,j,k
 	integer                   :: Nd_start(4)=0
 	integer                   :: Nd_count(4)=0
 	integer                   :: ii,jj,kk,offset
 	real(kind=8),allocatable  :: scratch(:,:,:)
 	real(kind=8)              :: h,deriv
     
    integer      :: locnx,locny,locnz
 	integer      :: fid,ncid,rcode  
 	integer      :: tVarID,s1VarID,s2VarID   
 	integer      :: divustarVarID,phiVarID,pdVarID
 	integer      :: uVarID,vVarID,wVarID
 	integer      :: ustarVarID,vstarVarID,wstarVarID
 	integer      :: uderivVarID,vderivVarID,wderivVarID
	integer      :: s1derivVarID,s2derivVarID
 	integer      :: do_bar,deriv_inc

 	if( .NOT. do_write(fid) ) then  !! no data to write
  		goto 999
 	endif
 	
 	
	locnx = array_size(JDIM,YBLOCK,myid)
	locny = array_size(IDIM,YBLOCK,myid)
	locnz = array_size(KDIM,YBLOCK,myid)

	!--------------------------------------------------- 
	!  open up the existing, initialized netcdf file
	!--------------------------------------------------- 
  	rcode=NF_OPEN(fullname(fid),NF_WRITE,ncid)
  	if(rcode.ne.NF_NOERR) then
   		write(0,*) '... ERROR OPENING NETCDF FILE: write_netcdf ',fid,trim(fullname(fid))
   		write(0,*) '... myid, rcode ',myid, rcode
   		stop
  	endif

	!------------------------------------------------------ 
	!   extract the variable ids given the variable names
	!    only look for ids for variables to be written
	!------------------------------------------------------   
  	rcode=NF_INQ_VARID(ncid,'time',tVarID)
  	if (rcode.ne.NF_NOERR) print *,myid,   &
      'NetCDF ERROR INQ_VARID -> time: ',rcode

  	if( variable_key(4,fid) /= 0 ) then
   		rcode=NF_INQ_VARID(ncid,'s1',s1VarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
      		'NetCDF ERROR INQ_VARID -> s1: ',rcode
  	endif
      
  	if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
   		rcode=NF_INQ_VARID(ncid,'s2',s2VarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
   			'NetCDF ERROR INQ_VARID -> s2: ',rcode
  	endif
      
  	if( variable_key(1,fid) == 1 ) then
		rcode=NF_INQ_VARID(ncid,'u',uVarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
       	'NetCDF ERROR INQ_VARID -> u: ',rcode
  	endif
      
  	if( variable_key(2,fid) == 1 ) then
    	rcode=NF_INQ_VARID(ncid,'v',vVarID)
     	if (rcode.ne.NF_NOERR) print *,myid,   &
       	 'NetCDF ERROR INQ_VARID -> v: ',rcode
  	endif
      
  	if( variable_key(3,fid) == 1 ) then
    	rcode=NF_INQ_VARID(ncid,'w',wVarID)
    	if (rcode.ne.NF_NOERR) print *,myid,   &
        	'NetCDF ERROR INQ_VARID -> w: ',rcode
  	endif
      
  	if( variable_key(6,fid) == 1 ) then        
    	rcode=NF_INQ_VARID(ncid,'divustar',divustarVarID)
    	if (rcode.ne.NF_NOERR) print *,myid,   &
      	  'NetCDF ERROR INQ_VARID -> divustar: ',rcode
  	endif
      
  	if( variable_key(7,fid) == 1 ) then        
   		rcode=NF_INQ_VARID(ncid,'phi',phiVarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
      	 'NetCDF ERROR INQ_VARID -> phi: ',rcode
  	endif
      
  	if( variable_key(8,fid) == 1 ) then        
   		rcode=NF_INQ_VARID(ncid,'pd',pdVarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
   	    'NetCDF ERROR INQ_VARID -> pd: ',rcode
  	endif
  	
  	if( variable_key(9,fid) == 1 ) then
		rcode=NF_INQ_VARID(ncid,'ustar',ustarVarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
       	'NetCDF ERROR INQ_VARID -> ustar: ',rcode
  	endif
  	
  	if( variable_key(10,fid) == 1 ) then
		rcode=NF_INQ_VARID(ncid,'vstar',vstarVarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
       	'NetCDF ERROR INQ_VARID -> vstar: ',rcode
  	endif
  	
  	if( variable_key(11,fid) == 1 ) then
		rcode=NF_INQ_VARID(ncid,'wstar',wstarVarID)
   		if (rcode.ne.NF_NOERR) print *,myid,   &
       	'NetCDF ERROR INQ_VARID -> wstar: ',rcode
  	endif
  	
  	if( t_secs >= t_start_child .and. do_child_grid .and. filename_root(fid) .NE. 'XYZ_child_' ) then
  		if( filename_root(fid)=='east' .or. filename_root(fid)=='west' ) then
  			rcode=NF_INQ_VARID(ncid,'u_x',uderivVarID)
  			rcode=NF_INQ_VARID(ncid,'v_x',vderivVarID)
  			rcode=NF_INQ_VARID(ncid,'w_x',wderivVarID)
  			rcode=NF_INQ_VARID(ncid,'s1_x',s1derivVarID)
  			if(write_s2_child_grid) rcode=NF_INQ_VARID(ncid,'s2_x',s2derivVarID) 
  			h = x(2)-x(1)
  		elseif( filename_root(fid)=='south' .or. filename_root(fid)=='north' ) then
  			rcode=NF_INQ_VARID(ncid,'u_y',uderivVarID)
  			rcode=NF_INQ_VARID(ncid,'v_y',vderivVarID)
  			rcode=NF_INQ_VARID(ncid,'w_y',wderivVarID)
  			rcode=NF_INQ_VARID(ncid,'s1_y',s1derivVarID)
  			if(write_s2_child_grid) rcode=NF_INQ_VARID(ncid,'s2_y',s2derivVarID)
  			h = y(2)-y(1)
  		elseif( filename_root(fid)=='bottom' .or. filename_root(fid)=='top' ) then
  			rcode=NF_INQ_VARID(ncid,'u_z',uderivVarID)
  			rcode=NF_INQ_VARID(ncid,'v_z',vderivVarID)
  			rcode=NF_INQ_VARID(ncid,'w_z',wderivVarID)
  			rcode=NF_INQ_VARID(ncid,'s1_z',s1derivVarID)
  			if(write_s2_child_grid) rcode=NF_INQ_VARID(ncid,'s2_z',s2derivVarID)
  			h = z(2)-z(1)
  		endif 		
  	endif
  	
  	
        
	!------------------------------------------------------------  
	!  set up count and start arrays for this file
	!------------------------------------------------------------
  	Nd_start(dimid_x(fid))=1
  	Nd_count(dimid_x(fid))=count(dimid_x(fid),fid)
       
  	Nd_start(dimid_y(fid))=1
 	Nd_count(dimid_y(fid))=count(dimid_y(fid),fid)
       
 	Nd_start(dimid_z(fid))=1
  	Nd_count(dimid_z(fid))=count(dimid_z(fid),fid)
       
  	Nd_start(dimid_t(fid))=time_counter(fid)
  	Nd_count(dimid_t(fid))=1   !! always 1 time slice written per call
  	offset = global_z_indices(START,YBLOCK,myid)-1
  
	!------------------------------------------------------------  
	!  write time to file
	!------------------------------------------------------------  
  	rcode=NF_PUT_VARA_DOUBLE(ncid,tVarID,Nd_start(dimid_t(fid)),        &
                         Nd_count(dimid_t(fid)),t_secs)
  	if(rcode.ne.NF_NOERR) print *,myid,'NetCDF ERROR tVarID: ',rcode
  	if( trim(mode(fid)) == 'append' ) time_counter(fid)=time_counter(fid)+1
  

	!------------------------------------------------------------  
	!  allocate a contiguous array for selected, 
	!  possibly strided output data
	!------------------------------------------------------------ 
  	allocate( scratch( Nd_count(dimid_x(fid)), &
         	           Nd_count(dimid_y(fid)), &
              	       Nd_count(dimid_z(fid)) ), stat=rcode )
  	if(rcode /= 0 ) stop 'problem allocating scratch in write_netcdf'
  
  
 	!------------------------------------------------------------
 	! Scalar 1  in dimensional units.
 	!------------------------------------------------------------
  	if( variable_key(4,fid) /= 0 ) then
   		if( write_s1_bar(fid) ) then
    		do_bar=1
   		else
    		do_bar=0
   		endif
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)         
           			scratch(ii,jj,kk) = do_bar*s1_bar(k+offset,1) + s1(j,i,k)                    
         			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     	ii=ii+1
    	enddo

 		!---------------------------------------------------
 		!   Each processor writes its local data to 
 		!   the appropriate locations in the netcdf file.
 		!---------------------------------------------------
   		rcode=NF_PUT_VARA_DOUBLE(ncid,s1VarID,Nd_start,Nd_count,scratch)
   		if(rcode.ne.NF_NOERR) then
    		print *,myid,'NetCDF ERROR: PUT_VARA_REAL scalar (1)',rcode
    		err_msg=NF_STRERROR(rcode)
    		print*,myid, err_msg
   		endif
	endif
      
	! Scalar 2  in dimensional units.
  	if( do_second_scalar .and. variable_key(5,fid) /= 0 ) then
   		if( write_s2_bar(fid) ) then
    		do_bar=1
   		else
    		do_bar=0
   		endif
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)         
           			scratch(ii,jj,kk) = do_bar*s2_bar(k+offset,1) + s2(j,i,k)        
         			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

 		!---------------------------------------------------
 		!   Each processor writes its local data to 
 		!   the appropriate locations in the netcdf file.
 		!---------------------------------------------------
   		rcode=NF_PUT_VARA_DOUBLE(ncid,s2VarID,Nd_start,Nd_count,scratch)
   		if(rcode.ne.NF_NOERR) then
    		print *,myid,'NetCDF ERROR: PUT_VARA_REAL scalar (2)',rcode
    		err_msg=NF_STRERROR(rcode)
    		print*,myid, err_msg
   		endif
	endif
      
	! u velocity component in dimensional units.
  	if( variable_key(1,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = u(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!  Each processor writes its local data to 
		!  the appropriate locations in the netcdf file.
   		rcode=NF_PUT_VARA_DOUBLE(ncid,uVarID,Nd_start,Nd_count,scratch)
   		if(rcode.ne.NF_NOERR) print *,myid,   &
      		'NetCDF ERROR: PUT_VARA_DOUBLE u',rcode
   	endif
      
	!      v velocity component in dimensional units.
  	if( variable_key(2,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = v(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!   Each processor writes its local data to 
		!   the appropriate locations in the netcdf file.
   		rcode=NF_PUT_VARA_DOUBLE(ncid,vVarID,Nd_start,Nd_count,scratch)
    	if(rcode.ne.NF_NOERR) print *,myid,   &
       		'NetCDF ERROR: PUT_VARA_DOUBLE v',rcode
    endif
      
	! w velocity component in dimensional units.
  	if( variable_key(3,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = w(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!   Each processor writes its local data to 
		!   the appropriate locations in the netcdf file.
       
    	rcode=NF_PUT_VARA_DOUBLE(ncid,wVarID,Nd_start,Nd_count,scratch)
    	if(rcode.ne.NF_NOERR) print *,myid,   &
       	'NetCDF ERROR PUT_VARA_DOUBLE: w ',rcode
	endif

	! div ustar in dimensional units.
	if( variable_key(6,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = div_u(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!   Each processor writes its local data to 
		!   the appropriate locations in the netcdf file.
    	rcode=NF_PUT_VARA_DOUBLE(ncid,divustarVarID,Nd_start,Nd_count,scratch)
    	if(rcode.ne.NF_NOERR) print *,myid,   &
       	'NetCDF ERROR PUT_VARA_DOUBLE: div ustar ',rcode
	endif

	! pressure soln in dimensional units.
  	if( variable_key(7,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = phi(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!   Each processor writes its local data to 
		!   the appropriate locations in the netcdf file.
    	rcode=NF_PUT_VARA_DOUBLE(ncid,phiVarID,Nd_start,Nd_count,scratch)
    	if(rcode.ne.NF_NOERR) print *,myid,   &
       	'NetCDF ERROR PUT_VARA_DOUBLE: phi ',rcode
    endif
    
    ! ustar velocity component in dimensional units.
  	if( variable_key(9,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = ustar(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!  Each processor writes its local data to 
		!  the appropriate locations in the netcdf file.
   		rcode=NF_PUT_VARA_DOUBLE(ncid,ustarVarID,Nd_start,Nd_count,scratch)
   		if(rcode.ne.NF_NOERR) print *,myid,   &
      		'NetCDF ERROR: PUT_VARA_DOUBLE ustar',rcode
   	endif
   	
   	! vstar velocity component in dimensional units.
  	if( variable_key(10,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = vstar(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!  Each processor writes its local data to 
		!  the appropriate locations in the netcdf file.
   		rcode=NF_PUT_VARA_DOUBLE(ncid,vstarVarID,Nd_start,Nd_count,scratch)
   		if(rcode.ne.NF_NOERR) print *,myid,   &
      		'NetCDF ERROR: PUT_VARA_DOUBLE vstar',rcode
   	endif
   	
   	! wstar velocity component in dimensional units.
  	if( variable_key(11,fid)==1 ) then
   		ii=1 ;  
    	do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     		jj=1;
      		do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       			kk=1;
        		do k=my_z0(fid),my_z1(fid),my_zinc(fid)
          			scratch(ii,jj,kk) = wstar(j,i,k)
          			kk=kk+1
        		enddo
       			jj=jj+1
      		enddo
     		ii=ii+1
    	enddo

		!  Each processor writes its local data to 
		!  the appropriate locations in the netcdf file.
   		rcode=NF_PUT_VARA_DOUBLE(ncid,wstarVarID,Nd_start,Nd_count,scratch)
   		if(rcode.ne.NF_NOERR) print *,myid,   &
      		'NetCDF ERROR: PUT_VARA_DOUBLE wstar',rcode
   	endif
   	
   	
   	if( t_secs >= t_start_child .and. do_child_grid .and. filename_root(fid) .NE. 'XYZ_child_' ) then
   	
  		if( filename_root(fid)=='east' .or. filename_root(fid)=='west' ) then
   	  			
   			deriv_inc = 1                             ! forward difference
   			if( my_x0(fid) == locnx ) deriv_inc = -1  ! backward difference needed
   			
   			! u deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( u(j,i+deriv_inc,k) - u(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,uderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE u_x',rcode
      			
      		! v deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( v(j,i+deriv_inc,k) - v(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,vderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE v_x',rcode
      			
      		! w deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( w(j,i+deriv_inc,k) - w(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,wderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE w_x',rcode
      			
      		! s1 deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( s1(j,i+deriv_inc,k) - s1(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,s1derivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE s1_x',rcode
      		
      		if( write_s2_child_grid ) then	
      			! s2 deriv
   				ii=1 ;  
    			do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     				jj=1;
      				do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       					kk=1;
        				do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        					deriv = deriv_inc*( s2(j,i+deriv_inc,k) - s2(j,i,k) )/h
          					scratch(ii,jj,kk) = deriv
          					kk=kk+1
        				enddo
       					jj=jj+1
      				enddo
     				ii=ii+1
    			enddo
   				rcode=NF_PUT_VARA_DOUBLE(ncid,s2derivVarID,Nd_start,Nd_count,scratch)
   				if(rcode.ne.NF_NOERR) print *,myid,   &
      				'NetCDF ERROR: PUT_VARA_DOUBLE s2_x',rcode  	
   			endif
   	
   		elseif( filename_root(fid)=='south' .or. filename_root(fid)=='north' ) then
   	  			
   			deriv_inc = 1                             ! forward difference
   			if( my_y0(fid) == locny ) deriv_inc = -1  ! backward difference needed
   			
   			! u deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( u(j+deriv_inc,i,k) - u(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,uderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE u_y',rcode
      			
      		! v deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( v(j+deriv_inc,i,k) - v(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,vderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE v_y',rcode
      			
      		! w deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( w(j+deriv_inc,i,k) - w(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,wderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE w_y',rcode
   	
   	
   			! s1 deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( s1(j+deriv_inc,i,k) - s1(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,s1derivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE s1_y',rcode
   	
   			if( write_s2_child_grid ) then
   				! s2 deriv
   				ii=1 ;  
    			do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     				jj=1;
      				do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       					kk=1;
        				do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        					deriv = deriv_inc*( s2(j+deriv_inc,i,k) - s2(j,i,k) )/h
          					scratch(ii,jj,kk) = deriv
          					kk=kk+1
        				enddo
       					jj=jj+1
      				enddo
     				ii=ii+1
    			enddo
   				rcode=NF_PUT_VARA_DOUBLE(ncid,s2derivVarID,Nd_start,Nd_count,scratch)
   				if(rcode.ne.NF_NOERR) print *,myid,   &
      				'NetCDF ERROR: PUT_VARA_DOUBLE s2_y',rcode
   			endif
   	
   		elseif( filename_root(fid)=='bottom' .or. filename_root(fid)=='top' ) then
   	  			
   			deriv_inc = 1                             ! forward difference
   			if( my_z0(fid) == locnz ) deriv_inc = -1  ! backward difference needed
   			
   			! u deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( u(j,i,k+deriv_inc) - u(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,uderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE u_z',rcode
      			
      		! v deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( v(j,i,k+deriv_inc) - v(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,vderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE v_z',rcode
      			
      		! w deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( w(j,i,k+deriv_inc) - w(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,wderivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE w_z',rcode
      			
      		! s1 deriv
   			ii=1 ;  
    		do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     			jj=1;
      			do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       				kk=1;
        			do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        				deriv = deriv_inc*( s1(j,i,k+deriv_inc) - s1(j,i,k) )/h
          				scratch(ii,jj,kk) = deriv
          				kk=kk+1
        			enddo
       				jj=jj+1
      			enddo
     			ii=ii+1
    		enddo
   			rcode=NF_PUT_VARA_DOUBLE(ncid,s1derivVarID,Nd_start,Nd_count,scratch)
   			if(rcode.ne.NF_NOERR) print *,myid,   &
      			'NetCDF ERROR: PUT_VARA_DOUBLE s1_z',rcode
      		
      		if( write_s2_child_grid ) then	
      			! s2 deriv
   				ii=1 ;  
    			do i=my_x0(fid),my_x1(fid),my_xinc(fid)
     				jj=1;
      				do j=my_y0(fid),my_y1(fid),my_yinc(fid)
       					kk=1;
        				do k=my_z0(fid),my_z1(fid),my_zinc(fid)
        					deriv = deriv_inc*( s2(j,i,k+deriv_inc) - s2(j,i,k) )/h
          					scratch(ii,jj,kk) = deriv
          					kk=kk+1
        				enddo
       					jj=jj+1
      				enddo
     				ii=ii+1
    			enddo
   				rcode=NF_PUT_VARA_DOUBLE(ncid,s2derivVarID,Nd_start,Nd_count,scratch)
   				if(rcode.ne.NF_NOERR) print *,myid,   &
      				'NetCDF ERROR: PUT_VARA_DOUBLE s2_z',rcode
      			endif
      			
      	endif
    endif
      
	!   Each processor closes the netcdf file.
    rcode=NF_CLOSE(ncid)
    if(rcode.ne.NF_NOERR) print *,myid,   &
      'NetCDF ERROR: closing file in write_netcdf',rcode
    deallocate( scratch )

999	continue
    call mpi_barrier(comm,rcode)  
 return
end subroutine write_netcdf



subroutine process_dimensions(fid,dimstr)
	use io_params
	implicit none
	integer                     :: fid,idim,jdim,ndims
	character(len=3)            :: dimstr
 
 	!-------------------------------------------------------------
 	! count the number of spatially varying indices up from 1
 	! if not varying, count down from 4
 	! nspace is the number of spatially varying indices
 	! these values are stored in dimid_x,dimid_y,dimid_z,dimid_t
 	!
 	! e.g. all spatial coordinates have > 1 output point
 	!      dimid_x,dimid_y,dimid_z,dimid_t ==> 1,2,3,4
 	!      nspace = 3
 	!
	 ! eg.  say y has only 1 value to be output
 	!      dimid_x,dimid_z,dimid_t,dimid_y ==> 1,2,3,4
 	!      nspace = 2
 	!
 	! eg.  say both y and x have only 1 value to be output
 	!      dimid_z,dimid_t,dimid_x,dimid_y ==> 1,2,3,4
 	!      nspace = 1
 	!
 	! eg.  say all of x,y,z have only 1 value to be output
 	!      dimid_t,dimid_x,dimid_y,dimid_z ==> 1,2,3,4
 	!      nspace = 0
 	!-------------------------------------------------------------
  
 	idim=1  !! count up for spatially varying dimensions
 	jdim=4  !! count down for dimensions w/ indices held constant
 	nspace(fid)=0    
  
 	if( ilocs(2,fid).ne.ilocs(1,fid) ) then  ! x nonsingleton output coord
  		nspace(fid)=nspace(fid)+1
  		dimid_x(fid) = idim
  		idim = idim + 1
 	else
  		dimid_x(fid) = jdim
  		jdim = jdim - 1
 	endif
   
 	if( jlocs(2,fid).ne.jlocs(1,fid) ) then  ! y nonsingleton output coord
  		nspace(fid)=nspace(fid)+1
  		dimid_y(fid) = idim
  		idim = idim + 1
 	else
  		dimid_y(fid) = jdim
  		jdim = jdim - 1
 	endif
   
 	if( klocs(2,fid).ne.klocs(1,fid) ) then  ! z nonsingleton output coord
 		 nspace(fid)=nspace(fid)+1
  		dimid_z(fid) = idim
  		idim = idim + 1
 	else
  		dimid_z(fid) = jdim
  		jdim = jdim - 1
 	endif
 	dimid_t(fid) = idim
 	ndims=idim   !! time + counted number of space dims

 	!------------------------------------------------------------------
 	! given # of space dimensions, decide which directory to write to
 	!------------------------------------------------------------------
 	if( ndims-1==0 ) then
 		dimstr='TS/'
 	elseif( ndims-1==1 ) then
 		dimstr='1D/'
 	elseif( ndims-1==2 ) then
 		dimstr='2D/'
 	elseif( ndims-1==3 ) then
		dimstr='3D/'
	endif
 return
end subroutine process_dimensions

subroutine construct_filename(fid,dimstr,topdir,istep,filename)
 use io_params
 use mpi_params,           only: myid,numprocs
 use decomposition_params, only: YBLOCK,proc_row,proc_col
 implicit none
 integer             :: fid
 integer             :: istep
 character(len=3)    :: dimstr
 character(len=80)   :: topdir
 character(len=80)   :: filename
 character(len=6)    :: cnum
 
 character(len=3)    :: c_row_id
 character(len=3)    :: c_col_id
 character(len=7)    :: cid
 

 write(unit = cnum, fmt = 100) istep
 write(unit = c_row_id,  fmt = 101) proc_row(YBLOCK,myid)
 write(unit = c_col_id,  fmt = 101) proc_col(YBLOCK,myid)
 
 cid = c_row_id//'-'//c_col_id
 !write(unit = cid,  fmt = 101) myid


 if( trim(mode(fid)) == 'new' ) then
  if(numprocs .gt. 1) then
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'_'//cnum//'_'//cid//'.nc'
  else
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'_'//cnum//'.nc'
  endif
 elseif( trim(mode(fid)) == 'append' ) then
  if(numprocs .gt. 1) then
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'_'//cid//'.nc'
  else
   filename = trim(topdir)//dimstr//trim(filename_root(fid))//'.nc'
  endif
 endif
 
100 format(I6.6)
101 format(I3.3)
end subroutine construct_filename



subroutine set_local_indices_count(fid)
	use io_params
	use decomposition_params
	use mpi_params,             only: myid
	implicit none
	integer                        :: fid
	integer                        :: offset(2)
	integer,external ::  count_vals

	do_write(fid)=.TRUE.   ! start by assuming there is data to write

	!----------------------------------------------------------------
	! Compute starting and ending indices into local YBLOCK storage
	! store counts and increments as well
	!----------------------------------------------------------------
 
	!----------------------------------------------------------------
	! YBLOCK ==> y coord is always local to myid
	!----------------------------------------------------------------
	my_y0(fid)   = jlocs(1,fid)
	my_y1(fid)   = jlocs(2,fid)
	my_yinc(fid) = jlocs(3,fid)
	my_ny(fid)   = count_vals(jlocs(1,fid),jlocs(2,fid),jlocs(3,fid))
	count(dimid_y(fid),fid) = my_ny(fid)
 
 
 	!----------------------------------------------------------------
 	! YBLOCK ==> x coord is distributed
 	!----------------------------------------------------------------
  	offset(1) = global_x_indices(START,YBLOCK,myid)
  	offset(2) = global_x_indices(END,YBLOCK,myid)
 
 	!-------------------------------------------------------------
 	! if desired data has end index lower than myid's start, or
 	! if desired data has start index higher than myid's end
 	! ===> no data on myid to write for this file set
 	!-------------------------------------------------------------
	if( ilocs(1,fid) > offset(2) .or. ilocs(2,fid) < offset(1) ) then
		do_write(fid)=.FALSE.
		return
	endif
 
 	!-------------------------------------------------------------
 	! myid has some data to be written, determine local
 	! starting and ending indices and number of points
 	!-------------------------------------------------------------
	my_x0(fid)   = maxval( (/ilocs(1,fid), offset(1)/) ) - offset(1) + 1
	my_x1(fid)   = minval( (/ilocs(2,fid), offset(2)/) ) - offset(1) + 1 
	my_xinc(fid) = ilocs(3,fid)
	my_nx(fid)   = count_vals(my_x0(fid),my_x1(fid),my_xinc(fid))
	count(dimid_x(fid),fid) = my_nx(fid)
  
  
 	!----------------------------------------------------------------
 	! YBLOCK ==> z coord is distributed
 	!----------------------------------------------------------------
	offset(1) = global_z_indices(START,YBLOCK,myid)
	offset(2) = global_z_indices(END,YBLOCK,myid)
 
 	!-------------------------------------------------------------
	 ! if desired data has end index lower than myid's start, or
	 ! if desired data has start index higher than myid's end
	 ! ===> no data on myid to write for this file set
	 !-------------------------------------------------------------
	if( klocs(1,fid) > offset(2) .or. klocs(2,fid) < offset(1) ) then
		do_write(fid)=.FALSE.
		return
	endif
 
 	!-------------------------------------------------------------
 	! myid has some data to be written, determine local
 	! starting and ending indices and number of points
 	!-------------------------------------------------------------
	my_z0(fid)   = maxval( (/klocs(1,fid), offset(1)/) ) - offset(1) + 1
	my_z1(fid)   = minval( (/klocs(2,fid), offset(2)/) ) - offset(1) + 1 
	my_zinc(fid) = klocs(3,fid)
	my_nz(fid)   = count_vals(my_z0(fid),my_z1(fid),my_zinc(fid)) 
	count(dimid_z(fid),fid) = my_nz(fid)   
 return 
 end subroutine set_local_indices_count
 
integer function count_vals(i0,i1,inc) result(numvals)
	implicit none
	integer, intent(in)  :: i0,i1,inc
	if( i1<i0 ) then
		numvals=0
		return
	endif
	numvals = floor( float(i1-i0)/inc ) + 1
end function count_vals




