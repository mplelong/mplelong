subroutine explicit_rhs
	use mpi_params,       only: myid
	use methods_params,   only: do_second_scalar
	use timing
	implicit none
		
	t_start_time_step = mpi_wtime()   ! start the clock at beginning of time step
	
	call rhs_u
	call rhs_v
	call rhs_w
	call rhs_s1
	if( do_second_scalar) then
		call rhs_s2
	endif	
 return
end subroutine explicit_rhs

subroutine rhs_u
	use mpi_params,                      only: myid
	use dependent_variables,             only: u,v,w
	use dimensional_scales,              only: f0
	use intermediate_variables,          only: explicit_rhs
	use etc,                             only: MM0
	implicit none
	integer,parameter                       :: id=1
	
	!------------------------------------------------------------------------
	! nonlinear term...
	!  udotgradf computes grad f and leaves f_x,f_y,f_z in tmpY(:,:,:,1-3)
	!------------------------------------------------------------------------
	call udotgradf( u, v, w, u, explicit_rhs(1,1,1,id,MM0) )
	explicit_rhs(:,:,:,id,MM0) = -explicit_rhs(:,:,:,id,MM0)     ! - uvec dot grad u
	
	!-----------------
	! rotation...
	!-----------------
	if( f0 .ne. 0.d0 ) then
		explicit_rhs(:,:,:,id,MM0) = explicit_rhs(:,:,:,id,MM0) + f0*v(:,:,:)
	endif	
 return
end subroutine rhs_u


subroutine rhs_v
	use mpi_params,                      only: myid
	use dependent_variables,             only: u,v,w
	use dimensional_scales,              only: f0
	use intermediate_variables,          only: explicit_rhs
	use etc,                             only: MM0
	implicit none
	integer,parameter                       :: id=2
	
	!-------------------------------------------------------------------
	! nonlinear term...
	!  udotgradf computes grad f and leaves results in tmpY(:,:,:,1-3)
	!-------------------------------------------------------------------
	call udotgradf( u, v, w, v, explicit_rhs(1,1,1,id,MM0) )
	explicit_rhs(:,:,:,id,MM0) = -explicit_rhs(:,:,:,id,MM0)     ! - uvec dot grad v

	!-----------------
	! rotation...
	!-----------------
	if( f0 .ne. 0.d0 ) then
		explicit_rhs(:,:,:,id,MM0) = explicit_rhs(:,:,:,id,MM0) - f0*u(:,:,:)
	endif
end subroutine rhs_v


subroutine rhs_w
	use mpi_params,                      only: myid
	use dependent_variables,             only: u,v,w,s1,scalar_kind
	use dimensional_scales,              only: rho0,g
	use intermediate_variables,          only: explicit_rhs
	use etc,                             only: MM0
	implicit none
	integer,parameter                       :: id=3
	
	!-------------------------------------------------------------------
	! nonlinear term...
	!  udotgradf computes grad f and leaves results in tmpY(:,:,:,1-3)
	!-------------------------------------------------------------------
	call udotgradf( u, v, w, w, explicit_rhs(1,1,1,id,MM0) )
	explicit_rhs(:,:,:,id,MM0) = -explicit_rhs(:,:,:,id,MM0)     ! - uvec dot grad w

	!-----------------
	! buoyancy...
	!-----------------
	if( scalar_kind(1) =='r' ) then
	 explicit_rhs(:,:,:,id,MM0) = explicit_rhs(:,:,:,id,MM0) - (g/rho0)*s1(:,:,:)    
	endif
end subroutine rhs_w



subroutine rhs_s1
!------------------------------------------------------------------------
! Set the rhs for the s1 equation. Don't include kappa*d2/dz2 s1_bar(z)
! because this is time independent and is integrated exactly in the
! time stepping routine.
!------------------------------------------------------------------------
	use mpi_params,                      only: myid
	use dependent_variables,             only: u,v,w,s1,s1_bar
	use intermediate_variables,          only: explicit_rhs
	use etc,                             only: MM0
	use methods_params,                  only: ambient_profile
	use decomposition_params
	implicit none
	integer,parameter                       :: id=4
	integer                                 :: i,j,k1,k2
	logical                                 :: debug=.TRUE.
	
	k1 = global_z_indices(START,YBLOCK,myid)
	k2 = global_z_indices(  END,YBLOCK,myid)
	!-------------------------------------------------------------------
	! nonlinear term...
	!  udotgradf computes grad f and leaves results in tmpY(:,:,:,1-3)
	!-------------------------------------------------------------------
	call udotgradf( u, v, w, s1, explicit_rhs(1,1,1,id,MM0) )
	explicit_rhs(:,:,:,id,MM0) = -explicit_rhs(:,:,:,id,MM0)     ! - uvec dot grad s1

	!-----------------
	! ambient profile
	!-----------------
	if( ambient_profile(1) ) then
		do i=1,array_size(JDIM,YBLOCK,myid)
			do j=1,array_size(IDIM,YBLOCK,myid)
				explicit_rhs(j,i,:,id,MM0) = explicit_rhs(j,i,:,id,MM0) - w(j,i,:)*s1_bar(k1:k2,2)
			enddo
		enddo		 
	endif
end subroutine rhs_s1

subroutine rhs_s2
!------------------------------------------------------------------------
! Set the rhs for the s2 equation. Don't include kappa*d2/dz2 s2_bar(z)
! because this is time independent and is integrated exactly in the
! time stepping routine.
!------------------------------------------------------------------------
	use mpi_params,                      only: myid
	use dependent_variables,             only: u,v,w,s2,s2_bar
	use intermediate_variables,          only: explicit_rhs
	use etc,                             only: MM0
	use methods_params,                  only: ambient_profile
	use decomposition_params
	implicit none
	integer,parameter                       :: id=5
	integer                                 :: i,j,k1,k2
	logical                                 :: debug=.TRUE.
	
	k1 = global_z_indices(START,YBLOCK,myid)
	k2 = global_z_indices(  END,YBLOCK,myid)	
	!-------------------------------------------------------------------
	! nonlinear term...
	!  udotgradf computes grad f and leaves results in tmpY(:,:,:,1-3)
	!-------------------------------------------------------------------
	call udotgradf( u, v, w, s2, explicit_rhs(1,1,1,id,MM0) )
	explicit_rhs(:,:,:,id,MM0) = -explicit_rhs(:,:,:,id,MM0)     ! - uvec dot grad s2

	!-----------------
	! ambient profile
	!-----------------
	if( ambient_profile(1) ) then
		do i=1,array_size(JDIM,YBLOCK,myid)
			do j=1,array_size(IDIM,YBLOCK,myid)
				explicit_rhs(j,i,:,id,MM0) = explicit_rhs(j,i,:,id,MM0) - w(j,i,:)*s2_bar(k1:k2,2)
			enddo
		enddo		 
	endif
end subroutine rhs_s2

