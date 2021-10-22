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
	call mudotgradf( u, v, w, u, explicit_rhs(1,1,1,id,MM0), id ) ! - uvec dot grad u
	
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
	call mudotgradf( u, v, w, v, explicit_rhs(1,1,1,id,MM0), id ) ! - uvec dot grad v

	!-----------------
	! rotation...
	!-----------------
	if( f0 .ne. 0.d0 ) then
		explicit_rhs(:,:,:,id,MM0) = explicit_rhs(:,:,:,id,MM0) - f0*u(:,:,:)
	endif
end subroutine rhs_v


subroutine rhs_w
	use mpi_params,                      only: myid,comm,ierr
	use dependent_variables,             only: u,v,w,s1,scalar_kind
	use dimensional_scales,              only: rho0,g
	use intermediate_variables,          only: explicit_rhs
	use etc,                             only: MM0
	use independent_variables,           only: Lz,z_FSRL,s1_z_BC,nx,ny,nz
	use decomposition_params
	implicit none
	integer,parameter                       :: id=3
	integer                                 :: k,icount
	integer, save                           :: locnz
	real(kind=8)                            :: dz,gamma,local_mean
	real(kind=8), allocatable               :: z(:)
	real(kind=8), allocatable, save         :: z_taper(:)
	real(kind=8), save                      :: rho_mean     ! spatial avg at t=0
	real(kind=8), external                  :: myexp
	logical, save                           :: first_entry=.TRUE.
	include "mpif.h"
		
	if( first_entry ) then					
		if( z_FSRL .and. s1_z_BC=='HOMOGENEOUS_NEUMANN') then
			locnz = array_size(KDIM,YBLOCK,myid)
			allocate( z(locnz), z_taper(locnz) )
			call get_my_zvals(z,YBLOCK,myid)
			dz = z(2)-z(1)
			gamma = 2.*dz
			do k=1,locnz
				z_taper(k) = 1.d0 - myexp(-((z(k)-0.d0)/gamma)**2) - myexp(-((z(k)-Lz)/gamma)**2)
			enddo
			deallocate(z)
		endif
		
!		if( scalar_kind(1)=='r' ) then
!			! calculate local spatial mean
!			local_mean = SUM( s1(:,:,:) )
!			call MPI_ALLREDUCE(local_mean,rho_mean,icount,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
!			rho_mean = rho_mean/(nx*ny*nz)
!		else
!			rho_mean = 0.d0
!		endif	
			
		first_entry=.FALSE.
	endif
	
	!-------------------------------------------------------------------
	! nonlinear term...
	!  udotgradf computes grad f and leaves results in tmpY(:,:,:,1-3)
	!-------------------------------------------------------------------
	call mudotgradf( u, v, w, w, explicit_rhs(1,1,1,id,MM0), id ) ! - uvec dot grad w

	!-----------------
	! buoyancy...
	!-----------------
	if( scalar_kind(1) =='r' ) then
!		explicit_rhs(:,:,:,id,MM0) = explicit_rhs(:,:,:,id,MM0) - (g/rho0)*(s1(:,:,:) - rho_mean)   
                explicit_rhs(:,:,:,id,MM0) = explicit_rhs(:,:,:,id,MM0) - (g/rho0)*(s1(:,:,:))
	endif
	
	! If necessary, taper acceleration due to buoyancy near top/bottom
	if( z_FSRL .and. s1_z_BC=='HOMOGENEOUS_NEUMANN') then 
		do k=1,locnz
			explicit_rhs(:,:,k,id,MM0) = explicit_rhs(:,:,k,id,MM0)*z_taper(k)
		enddo  
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
	call mudotgradf( u, v, w, s1, explicit_rhs(1,1,1,id,MM0), id ) ! - uvec dot grad s1

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
	call mudotgradf( u, v, w, s2, explicit_rhs(1,1,1,id,MM0), id ) ! - uvec dot grad s2

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

