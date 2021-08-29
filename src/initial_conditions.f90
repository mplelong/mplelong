subroutine InitialConditions  
	use etc,                    only: logfile
	use independent_variables,  only: x,y,z,t0,tf
	use methods_params,         only: do_second_scalar,ambient_profile,restart
	use methods_params,         only: rs_basename,subtract_s1_bar,subtract_s2_bar,add_restart_time
	use decomposition_params 
	use intermediate_variables
	use dependent_variables
	use dimensional_scales
	use mpi_params
	implicit none 
 
 	include 'netcdf.inc'
	real(kind=8),allocatable   :: my_xvals(:)
	real(kind=8),allocatable   :: my_yvals(:)
	real(kind=8),allocatable   :: my_zvals(:)
	integer                    :: id,locnx,ny,locnz,i,j,k,kg
	
	integer                    ::  ncid,varid,i_proc,j_proc
	character(len=4)           ::  cid
	character(len=3)           ::  iproc,jproc
	character(len=80)          ::  ncfile
	integer                    ::  nc_start(4),count(4)
	logical                    ::  debug=.FALSE.
	real(kind=8)               ::  restart_time
	real(kind=8),allocatable   ::  tmp(:,:,:,:) 
 
  
	if(myid==0) then
		write(0,*) ' ................'
		write(0,*) ' ................     hello world from InitialConditions'
		open(1,file=logfile,position='append') 
		write(1,*) '  '
		write(1,*) '  '
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '                  InitialConditions Report:'
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '  '
		close(1)
	endif
	
	locnx = array_size(JDIM,YBLOCK,myid)   
	ny = array_size(IDIM,YBLOCK,myid)  
	locnz = array_size(KDIM,YBLOCK,myid)  

  	
  	if( .not. restart )	then
		!------------------------------------------------
		!    YBLOCK arranged as data(y,x,z)
		!    data(ny,locnx,locnz)
		!------------------------------------------------
		allocate( my_xvals(locnx) ) ; my_xvals(:)=0.d0
		allocate( my_yvals(ny) ) ; my_yvals(:)=0.d0
		allocate( my_zvals(locnz) ) ; my_zvals(:)=0.d0
	
		call get_my_xvals( my_xvals, YBLOCK, myid )
		call get_my_yvals( my_yvals, YBLOCK, myid )
		call get_my_zvals( my_zvals, YBLOCK, myid )
   		 
 		id=1
 		call user_ics(my_xvals,my_yvals,my_zvals,u,id,locnx,ny,locnz)
 
 		id=2
		call user_ics(my_xvals,my_yvals,my_zvals,v,id,locnx,ny,locnz)

		id=3
 		call user_ics(my_xvals,my_yvals,my_zvals,w,id,locnx,ny,locnz)
 
		id=4
		call user_ics(my_xvals,my_yvals,my_zvals,s1,id,locnx,ny,locnz) 
 
		if( do_second_scalar ) then 
  			id=5
  			call user_ics(my_xvals,my_yvals,my_zvals,s2,id,locnx,ny,locnz)
 		endif
 
		deallocate( my_xvals )
		deallocate( my_yvals )
		deallocate( my_zvals ) 
 
 	elseif( restart ) then
 	
 		!------------------------------------------------------
 		!   construct the filename for this processor to read
 		!------------------------------------------------------
  		write(unit=iproc,fmt='(I3.3)') proc_row(YBLOCK,myid)    ! 3 digit char string for iproc
  		write(unit=jproc,fmt='(I3.3)') proc_col(YBLOCK,myid)    ! 3 digit char string for jproc
		ncfile=trim(rs_basename)//iproc//'-'//jproc//'.nc'      ! i.e. RESTART/XYZ_123607_000-071.nc 
		
		!-----------------------------------------------------------------------------------
		! For renc_start files created by flow_solve itself OR create_renc_start_files.py
		! ncdump -h reports e.g.  v(timedimension,kdimension, jdimension, idimension) in 3d
		! but the order of the indices needs to be reversed here for fortran ==>(x,y,z,t)
		!-----------------------------------------------------------------------------------
		nc_start=(/1,1,1,1/)
		count=(/locnx,ny,locnz,1/)          ! ny (global) locnx, locnz for myid
		allocate( tmp(locnx,ny,locnz,1) )   ! tmp array used for reading in the data (NB not YBLOCK)

		!------------------------------------------------------
  		!  open the file, check for error
  		!------------------------------------------------------
		ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
		if(ierr.ne.NF_NOERR) then
			write(0,*) '... ERROR OPENING NETCDF FILE: user_ics ',trim(ncfile)
			write(0,*) '... myid, ierr ',myid, ierr
			stop
		endif
		
		ierr=NF_INQ_VARID(ncid,'time',varid)
		ierr = NF_GET_VARA_DOUBLE(ncid,varid,1,1,restart_time)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> restart_time: ',ierr
			stop
		endif
		
		if( add_restart_time ) then    ! add the time to the user vals in problem_params
			t0 = t0 + restart_time
			tf = tf + restart_time
		endif
				

		!------------------------------------------------------
		!  extract variable id, check for error
		!------------------------------------------------------
		ierr=NF_INQ_VARID(ncid,'u',varid)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR INQ_VARID -> u: ',ierr
			stop
		endif

		!------------------------------------------------------
		!  read the corresponding variable (dimensional)
		!------------------------------------------------------
		ierr = NF_GET_VARA_DOUBLE(ncid,varid,nc_start,count,tmp)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,' NetCDF ERROR NF_GET_VARA_REAL -> u: ',ierr
			write(0,*) myid,trim(ncfile),nc_start,count
			stop
		endif

	   !----------------------------------------
	   !  u
	   !----------------------------------------
		do k=1,locnz
			do i=1,locnx
				do j=1,ny
					u(j,i,k) = tmp(i,j,k,1)   ! 3d XYZ file
				enddo
			enddo
		enddo
		if( myid==0 .and. debug ) write(0,*) 'user_ics: u vals read and stored'
   

		!------------------------------------------------------
		!  extract variable id, check for error
		!------------------------------------------------------
		ierr=NF_INQ_VARID(ncid,'v',varid)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR INQ_VARID -> v: ',ierr
			stop
		endif

		!------------------------------------------------------
		!  read the corresponding variable (dimensional)
		!------------------------------------------------------
		ierr = NF_GET_VARA_DOUBLE(ncid,varid,nc_start,count,tmp)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> v: ',ierr
			stop
		endif

	   !----------------------------------------
	   !  v
	   !----------------------------------------
		do k=1,locnz
			do i=1,locnx
				do j=1,ny
					v(j,i,k) = tmp(i,j,k,1)   ! 3d XYZ file
				enddo
			enddo
		enddo
		if( myid==0 .and. debug ) write(0,*) 'user_ics: v vals read and stored'
   

		!------------------------------------------------------
		!  extract variable id, check for error
		!------------------------------------------------------
		ierr=NF_INQ_VARID(ncid,'w',varid)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR INQ_VARID -> v: ',ierr
			stop
		endif

		!------------------------------------------------------
		!  read the corresponding variable (dimensional)
		!------------------------------------------------------
		ierr = NF_GET_VARA_DOUBLE(ncid,varid,nc_start,count,tmp)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> w: ',ierr
			stop
		endif

	   !----------------------------------------
	   !  w
	   !----------------------------------------
		do k=1,locnz
			do i=1,locnx
				do j=1,ny
					w(j,i,k) = tmp(i,j,k,1)   ! 3d XYZ file
				enddo
			enddo
		enddo
		if( myid==0 .and. debug ) write(0,*) 'user_ics: w vals read and stored'
   
		!------------------------------------------------------
		!  extract variable id, check for error
		!------------------------------------------------------
		ierr=NF_INQ_VARID(ncid,'s1',varid)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR INQ_VARID -> s1: ',ierr
			stop
		endif

		!------------------------------------------------------
		!  read the corresponding variable (dimensional)
		!------------------------------------------------------
		ierr = NF_GET_VARA_DOUBLE(ncid,varid,nc_start,count,tmp)
		if (ierr.ne.NF_NOERR) then
			write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> s1: ',ierr
			stop
		endif

	   !----------------------------------------
	   !  s1
	   !----------------------------------------
		do k=1,locnz
			do i=1,locnx
				do j=1,ny
					s1(j,i,k) = tmp(i,j,k,1)   ! 3d XYZ file
				enddo
			enddo
		enddo
		if( myid==0 .and. debug ) write(0,*) 'user_ics: s1 vals read and stored'
 
   		!------------------------------------------------------
   		!  s1 + s1_bar in netcdf file ==> subtract s1_bar
   		!   whether to do this depends on how the
   		!   original output was done, which is optional
   		!------------------------------------------------------
   		if( subtract_s1_bar ) then
    		do k=1,locnz
    			kg = global_z_indices(START,YBLOCK,myid) + k - 1  ! global k index
     			s1(:,:,k) = s1(:,:,k) - s1_bar(kg,1)
    		enddo
   		endif
   		if( myid==0 .and. debug ) write(0,*) 'user_ics: s1 vals read and stored'
   

 		if( do_second_scalar ) then
			!------------------------------------------------------
			!  extract variable id, check for error
			!------------------------------------------------------
			ierr=NF_INQ_VARID(ncid,'s2',varid)
			if (ierr.ne.NF_NOERR) then
				write(0,*) myid,'NetCDF ERROR INQ_VARID -> s2: ',ierr
				stop
			endif

			!------------------------------------------------------
			!  read the corresponding variable (dimensional)
			!------------------------------------------------------
			ierr = NF_GET_VARA_DOUBLE(ncid,varid,nc_start,count,tmp)
			if (ierr.ne.NF_NOERR) then
				write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> s2: ',ierr
				stop
			endif

	   		!----------------------------------------
	   		!  s2
	   		!----------------------------------------
			do k=1,locnz
				do i=1,locnx
					do j=1,ny
						s2(j,i,k) = tmp(i,j,k,1)   ! 3d XYZ file
					enddo
				enddo
			enddo
			if( myid==0 .and. debug ) write(0,*) 'user_ics: s2 vals read and stored'
 
   			!------------------------------------------------------
   			!  s2 + s2_bar in netcdf file ==> subtract s2_bar
   			!   whether to do this depends on how the
   			!   original output was done, which is optional
   			!------------------------------------------------------
   			if( subtract_s2_bar ) then
    			do k=1,locnz
    				kg = global_z_indices(START,YBLOCK,myid) + k - 1  ! global k index
     				s2(:,:,k) = s2(:,:,k) - s2_bar(kg,1)
    			enddo
   			endif
   			if( myid==0 .and. debug ) write(0,*) 'user_ics: s2 vals read and stored'
		endif
		
 		deallocate( tmp )
 		ierr = NF_CLOSE(ncid)
 	
 	endif  ! end restart logic block
 
 	if(myid==0) then
  		open(1,file=logfile,position='append')
  			if( do_second_scalar ) then
  				write(1,*) ' ................      initialization of YBLOCK arrays u,v,w,s1,s2:  '
  			else
  				write(1,*) ' ................      initialization of YBLOCK arrays u,v,w,s1:     '
  			endif
  			write(1,*) ' ................      initialization via flow_solve_user_routines.f90/user_ics '
  			write(1,*) ' -----> InitialConditions routine exiting normally  <---------- '
  		close(1)
	endif
 
	call mpi_barrier(comm,ierr)
 return
end subroutine InitialConditions
