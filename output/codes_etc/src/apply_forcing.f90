subroutine apply_forcing
!---------------------------------------------------------------
!  Call the user routine user_forcing, passing in local
!  x, y and z values. Add user's forcing to rhs(:,:,:,id,MM0).
!---------------------------------------------------------------
	use mpi_params,                          only: myid
	use decomposition_params
	use independent_variables,               only: x,y,z,t_secs
	use intermediate_variables,              only: tmpY,explicit_rhs
	use methods_params,                      only: do_second_scalar,forcing_key,do_forcing
	use etc,                                 only: MM0 
	implicit none
	integer                                      :: id,k 
	real(kind=8),allocatable,save                :: xvals(:),yvals(:),zvals(:)
	integer,save                                 :: npts(3),npvs
	logical,save                                 :: first_entry=.TRUE.
  
	if( do_forcing ) then
 
		if( first_entry ) then
			npts(1) = array_size(JDIM,YBLOCK,myid)  ! # of x pts in YBLOCK
			npts(2) = array_size(IDIM,YBLOCK,myid)  ! # of y pts in YBLOCK
			npts(3) = array_size(KDIM,YBLOCK,myid)  ! # of z pts in YBLOCK
			if( do_second_scalar ) then
   				npvs=5
  			else
   				npvs=4
  			endif
  
			allocate( xvals( npts(1) ) )
			allocate( yvals( npts(2) ) )
			allocate( zvals( npts(3) ) )
		
			call get_my_xvals( xvals, YBLOCK, myid )
			call get_my_yvals( yvals, YBLOCK, myid )
			call get_my_zvals( zvals, YBLOCK, myid )
		
			first_entry=.FALSE.
		endif
 
  
 		!----------------------------------------- 
 		! First pass, assume each field is forced 
 		! but then let user have finer control over
 		! which variables need to be forced.
 		! Time independent forcing fields are not
 		! saved, calls are needed each time step.
 		!-----------------------------------------
		do id=1,npvs
			if( forcing_key(id) ) then 
				call user_forcing(xvals,yvals,zvals,t_secs,tmpY(:,:,:,1),id,    &                       
                   	              npts(1),npts(2),npts(3),forcing_key(id)) 
				!------------------------------------
				! add to previously computed values
				!------------------------------------
				explicit_rhs(:,:,:,id,MM0) = explicit_rhs(:,:,:,id,MM0) + tmpY(:,:,:,1) 
  			endif
 		enddo  
 
	endif  ! end if do_forcing block
   
 return
end subroutine apply_forcing
