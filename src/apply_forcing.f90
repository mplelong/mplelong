subroutine apply_forcing
!---------------------------------------------------------------
!  Call the user routine user_forcing, passing in local
!  x, y and z values. Add user's forcing to rhs(:,:,:,id,MM0).
!---------------------------------------------------------------
	use mpi_params,                          only: myid
	use decomposition_params
	use independent_variables,               only: x,y,z,t_secs,z_FSRL,s1_z_BC,Lz
	use intermediate_variables,              only: tmpY,explicit_rhs
	use methods_params,                      only: do_second_scalar,forcing_key,do_forcing
	use etc,                                 only: MM0 
	implicit none
	integer                                      :: id,k 
	real(kind=8),allocatable,save                :: xvals(:),yvals(:),zvals(:),z_taper(:)
	integer,save                                 :: npts(3),npvs
	real(kind=8)                                 :: dz,gamma
	real(kind=8), external                       :: myexp
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
			
			if( z_FSRL ) then
				allocate( z_taper(npts(3)) )
				dz = zvals(2)-zvals(1)
				gamma = 2.d0*dz
				do k=1,npts(3)
					z_taper(k) = 1.d0 - myexp(-((zvals(k)-0.d0)/gamma)**2) - myexp(-((zvals(k)-Lz)/gamma)**2)
				enddo
			endif
		
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
				
				!---------------------------------------------------------------------
 				! If necessary, taper vertical acceleration near top/bottom
 				!  NB nontraditional Coriolis term also needs tapering independent
 				! of the buoyancy term
 				!---------------------------------------------------------------------
				if( id==3 .and. z_FSRL .and. s1_z_BC=='HOMOGENEOUS_NEUMANN') then 
					do k=1,npts(3)
						explicit_rhs(:,:,k,id,MM0) = explicit_rhs(:,:,k,id,MM0)*z_taper(k)
					enddo  
				endif  
  			endif
 		enddo 
 		
 		
 
	endif  ! end if do_forcing block
   
 return
end subroutine apply_forcing
