subroutine pressure_projection 
!------------------------------------------------------------
!   main routine controlling the  pressure projection
!------------------------------------------------------------
	use mpi_params,               only: myid
	use intermediate_variables,   only: div_u,phi,ustar,vstar,wstar
	use independent_variables,    only: Lx,Ly,Lz,x_periodic,y_periodic,z_periodic
	use etc,                      only: istep,istart 
	use decomposition_params
	implicit none
	real(kind=8),allocatable,save    :: x(:),y(:),z(:)   ! local values
	integer, save                    :: locnx,ny,locnz,npts(3)
	integer                          :: i,j,k
	character(len=80)                :: dir
	logical, save                    :: first_entry=.TRUE., test_flag=.FALSE.
	
	if( first_entry) then
		npts(:) = 4                                ! near-boundary smoothing scale		
		locnx = array_size(JDIM,YBLOCK,myid)
		ny    = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)		
		!--------------------
		! local x,y,z vals
		!--------------------
		allocate( x(locnx), y(ny), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)		
		first_entry = .FALSE.
	endif
	
  
	!---------------------------------------------------------------
	! Update u* with boundary information. This enables homogeneous
	! boundary conditions to be used for the pressure variable phi.
	!---------------------------------------------------------------
	call update_ustar
 
	!-------------------------------------------------------------
	! Compute divergence (u*), uses tmpY(:,:,:,6) for work space.
	!-------------------------------------------------------------
	call divergence(ustar,vstar,wstar,div_u)
	
	dir='xy'
	call boundary_smooth(div_u,dir,npts)     ! this helps with w at east and west boundaries
	
	if( istep==istart ) then
	
		if(myid==0) then
			i=1;k=locnz/2
			open(1,file='output/debug_data/div_ustar_of_y')
				do j=1,ny
					write(1,*) y(j),div_u(j,i,k)
				enddo
			close(1)
		endif
		
		if( z(1)==0.d0 ) then
			i=1;j=ny/2 + 1
 			open(1,file='output/debug_data/div_ustar_of_z_B')
 				do k=1,locnz
 					write(1,*) z(k),div_u(j,i,k)
 				enddo
 			close(1)
 		endif
 		 		
 		if( z(locnz)==Lz ) then
 			i=1;j=ny/2 + 1
 			open(1,file='output/debug_data/div_ustar_of_z_T')
 				do k=1,locnz
 					write(1,*) z(k),div_u(j,i,k)
 				enddo
 			close(1)
 		endif
 		
	endif
     
	!-------------------------------------------------------------
	! Solve the poisson eqn using cos series expansion. 
	!-------------------------------------------------------------
	call poisson_solver(test_flag)
  
	!-------------------------------------------------------------
	! Compute the gradient of phi and project.
	!-------------------------------------------------------------
	call project
   
 return 
end subroutine pressure_projection


subroutine project
	!-------------------------------------------------------------
	! on entry, u*,v* and w* are stored in u, v and w arrays
	!   u <- u* -dt*GradPhi(1)  
	!   v <- v* -dt*GradPhi(2)
	!   w <- w* -dt*GradPhi(3)
	!-------------------------------------------------------------
	use mpi_params,                      only: myid
	use dependent_variables,             only: u,v,w
	use independent_variables,           only: Lx,Ly,Lz
	use intermediate_variables,          only: phi,tmpY,ustar,vstar,wstar
	use independent_variables,           only: dt
	use etc,                             only: istep,istart
	use decomposition_params
	implicit none 
	integer,parameter                            :: Q0=0
	integer                                      :: i, j, k
	integer, save                                :: locnx, ny, locnz
	real(kind=8),allocatable,dimension(:),save   :: x, y, z   ! local portions
	logical,save                                 :: first_entry=.TRUE.
 
	if( first_entry ) then	
		locnx = array_size(JDIM,YBLOCK,myid)
		ny    = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)		
		!--------------------
		! local x,y,z vals
		!--------------------
		allocate( x(locnx), y(ny), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)
		first_entry=.FALSE.
	endif
 
	!--------------------------------------------------------------
	! Compute grad phi...by design phi satisfies homogeneous
	! Neumann conditions so term by term differentiation of
	! it's cos series is ok. Use Q=0 to suppress Bernoulli scheme.
	!-------------------------------------------------------------- 
	call gradient(phi,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),Q0)
 
	!------------------------------------------------------------
	! update the velocity components
	!------------------------------------------------------------
	u = ustar - dt * tmpY(:,:,:,1)      ! u <== u* - dt phi_x
	v = vstar - dt * tmpY(:,:,:,2)      ! v <== v* - dt phi_y
	w = wstar - dt * tmpY(:,:,:,3)      ! w <== w* - dt phi_z
	
	if( istep==istart ) then
	
		if(myid==0) then
			i=1 ; k=locnz/2
 			open(1,file='output/debug_data/grad_phi_of_y')
 				do j=1,ny
 					write(1,*) y(j),tmpY(j,i,k,1),tmpY(j,i,k,2),tmpY(j,i,k,3)
 				enddo	
 			close(1)
 		endif		
 		if( z(1)==0.d0 ) then
 			open(1,file='output/debug_data/grad_phi_of_z_B')
 				i = 1; j = ny/2 + 1
 				do k=1,locnz
 					write(1,*) z(k),tmpY(j,i,k,1),tmpY(j,i,k,2),tmpY(j,i,k,3)
 				enddo
 			close(1)
 		endif 		 		
 		if( z(locnz)==Lz ) then
 			i = 1; j = ny/2 + 1
 			open(1,file='output/debug_data/grad_phi_of_z_T')
 				do k=1,locnz
 					write(1,*) z(k),tmpY(j,i,k,1),tmpY(j,i,k,2),tmpY(j,i,k,3)
 				enddo
 			close(1)
 		endif
 		
 		if(myid==0) then
			i=1 ; k=locnz/2
 			open(1,file='output/debug_data/u_projected_of_y')
 				do j=1,ny
 					write(1,*) y(j),u(j,i,k),v(j,i,k),w(j,i,k)
 				enddo	
 			close(1)
 		endif		
 		if( z(1)==0.d0 ) then
 			open(1,file='output/debug_data/u_projected_of_z_B')
 				i = 1; j = ny/2 + 1
 				do k=1,locnz
 					write(1,*) z(k),u(j,i,k),v(j,i,k),w(j,i,k)
 				enddo
 			close(1)
 		endif 		 		
 		if( z(locnz)==Lz ) then
 			i = 1; j = ny/2 + 1
 			open(1,file='output/debug_data/u_projected_of_z_T')
 				do k=1,locnz
 					write(1,*) z(k),u(j,i,k),v(j,i,k),w(j,i,k)
 				enddo
 			close(1)
 		endif
 		
	endif
  
 return  
end subroutine project


subroutine update_ustar
	!-------------------------------------------------------------
	! on entry, u*,v* and w* are stored in u, v and w arrays
	! add boundary info to u*, v* and w* arrays
	!-------------------------------------------------------------
	use mpi_params,                           only: myid
	use decomposition_params
	use dependent_variables,                  only: u,v,w,s1,s2
	use independent_variables,                only: Lx,Ly,Lz,x_periodic,y_periodic,z_periodic
	use intermediate_variables,               only: ustar,vstar,wstar,tmpY
	use boundary_data
	use etc,                                  only: istep,istart
	implicit none
 
 	integer, save                                :: locnx, ny, locnz, p, npts(3)
 	integer                                      :: id, i, j, k
 	real(kind=8)                                 :: gamma, a, b
	real(kind=8),allocatable,dimension(:),save   :: x, y, z   ! local portions
	real(kind=8),external                        :: myexp     ! to avoid large negative args to intrinsic exp()
	logical,save                                 :: first_entry=.TRUE.
 
	if( first_entry ) then	
		p = 4       ! exponent
		npts(:) =  (/3,3,6/)   ! gamma = npts*h, h=step size   2,2,6 is also good
		!--------------------
		! array map
		!--------------------
		! psi_u => tmpY(:,:,:,4)
		! psi_v => tmpY(:,:,:,5)
		! psi_w => tmpY(:,:,:,6)
		
		locnx = array_size(JDIM,YBLOCK,myid)
		ny    = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		
		!--------------------
		! local x,y,z vals
		!--------------------
		allocate( x(locnx), y(ny), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)
		
		!---------------------------------------------------
		! allocate near-boundary approximate step functions
		!---------------------------------------------------
		allocate( step_e(locnx), step_w(locnx) )
		allocate( step_s(ny)   , step_n(ny)    )
		allocate( step_b(locnz), step_t(locnz) )
		step_e=0.d0; step_w=0.d0; step_s=0.d0; step_n=0.d0; step_b=0.d0; step_t=0.d0
		
		!---------------------------------------------------
		! define near-boundary approximate step functions
		!---------------------------------------------------
		if( locnx > 1 ) then
			a = 0.d0
			b = Lx
			gamma = npts(1)*( x(2)-x(1) )
			do i=1,locnx
				step_e(i) = myexp(-((x(i)-a)/gamma)**p)   ! approx step function near east bdry
				step_w(i) = myexp(-((x(i)-b)/gamma)**p)   ! approx step function near west bdry
			enddo
		endif
		
		if( ny > 1 ) then
			a = 0.d0
			b = Ly
			gamma = npts(2)*( y(2)-y(1) )
			do j=1,ny
				step_s(j) = myexp(-((y(j)-a)/gamma)**p)   ! approx step function near south bdry
				step_n(j) = myexp(-((y(j)-b)/gamma)**p)   ! approx step function near north bdry
			enddo
		endif
		
		if( locnz > 1 ) then
			a = 0.d0
			b = Lz
			gamma = npts(3)*( z(2)-z(1) )
			do k=1,locnz
				step_b(k) = myexp(-((z(k)-a)/gamma)**p)   ! approx step function near bottom bdry
				step_t(k) = myexp(-((z(k)-b)/gamma)**p)   ! approx step function near top bdry
			enddo
		endif
		first_entry=.FALSE.
	endif
	
	!------------------------------------------------------------------------------	
	! construct an array that matches the required Dirichlet bcs at east/west
	!------------------------------------------------------------------------------
	id = 1  ! u
	if( locnx > 1 .and. .NOT. x_periodic ) then
		tmpY(:,:,:,4) = 0.d0
		if( x(1) == 0. ) then     ! myid owns east boundary, update psi_u
			do k=1,locnz
				do j=1,ny
					a = east_vals(j,k,id) - ustar(j,1,k)
					tmpY(j,:,k,4) =  a*step_e(:)
				enddo
			enddo
		endif
		if( x(locnx) == Lx ) then     ! myid owns west boundary, update psi_u
			do k=1,locnz
				do j=1,ny
					b = west_vals(j,k,id) - ustar(j,locnx,k)
					tmpY(j,:,k,4) = tmpY(j,:,k,4) + b*step_w(:)
				enddo
			enddo
		endif
	endif
	
	!------------------------------------------------------------------------------	
	! construct an array that matches the required Dirichlet bcs at south/north
	!------------------------------------------------------------------------------
	id = 2  ! v  YBLOCK, myid owns both y=0 and y=Ly
	if( ny > 1 .and. .NOT. y_periodic ) then
		tmpY(:,:,:,5) = 0.d0
		do k=1,locnz
			do i=1,locnx
				a = south_vals(i,k,id) - vstar(1,i,k)
				b = north_vals(i,k,id) - vstar(ny,i,k)
				tmpY(1:ny,i,k,5) = a*step_s(1:ny) + b*step_n(1:ny)
			enddo
		enddo
	endif
	
	!------------------------------------------------------------------------------	
	! construct an array that matches the required Dirichlet bcs at bottom/top
	!------------------------------------------------------------------------------
	id = 3  ! w
	if( locnz > 1 .and. .NOT. z_periodic ) then
		tmpY(:,:,:,6) = 0.d0
		if( z(1) == 0. ) then     ! myid owns bottom boundary, update psi_w
			do i=1,locnx
				do j=1,ny
					a = bottom_vals(i,j,id) - wstar(j,i,1)
					tmpY(j,i,:,6) =  a*step_b(:)
				enddo
			enddo
		endif
		if( z(locnz) == Lz ) then     ! myid owns top boundary, update psi_w
			do i=1,locnx
				do j=1,ny
					b = top_vals(i,j,id) - wstar(j,i,locnz)
					tmpY(j,i,:,6) = tmpY(j,i,:,6) + b*step_t(:)
				enddo
			enddo
		endif
	endif
	
	!---------------------------------------------
	! write some verification data at 1st step
	! x,y,z are local values; locnx,ny,locnz
	!---------------------------------------------
	if( istep==istart ) then
		if(myid==0 .and. ny > 1) then
			i=1 ; k=locnz/2
 			open(1,file='output/debug_data/psi_v_of_y')
 				do j=1,ny
 					write(1,*) y(j),tmpY(j,i,k,5)
 				enddo	
 			close(1)
 		endif
 		
 		if(myid==0 .and. locnx > 1) then
			j=1 ; k=locnz/2
 			open(1,file='output/debug_data/psi_u_of_x')
 				do i=1,locnx
 					write(1,*) x(i),tmpY(j,i,k,4)
 				enddo	
 			close(1)
 		endif
 	
 		
 		if( z(1)==0.d0 ) then
 			open(1,file='output/debug_data/psi_w_of_z_B')
 				i = 1; j = ny/2 + 1
 				do k=1,locnz
 					write(1,*) z(k),tmpY(j,i,k,6)
 				enddo
 			close(1)
 		endif
 		 		
 		if( z(locnz)==Lz ) then
 			i = 1; j = ny/2 + 1
 			open(1,file='output/debug_data/psi_w_of_z_T')
 				do k=1,locnz
 					write(1,*) z(k),tmpY(j,i,k,6)
 				enddo
 			close(1)
 		endif
 		
 		if(myid==0 .and. ny > 1) then
			i=1 ; k=locnz/2
 			open(1,file='output/debug_data/ustar_of_y')
 				do j=1,ny
 					write(1,*) y(j),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo	
 			close(1)
 		endif
 		
 		if(myid==0 .and. locnx > 1) then
			j=1 ; k=locnz/2
 			open(1,file='output/debug_data/ustar_of_x')
 				do i=1,locnx
 					write(1,*) x(i),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo	
 			close(1)
 		endif
 	
 		
 		if( z(1)==0.d0 ) then
 			open(1,file='output/debug_data/ustar_of_z_B')
 				i = 1; j = ny/2 + 1
 				do k=1,locnz
 					write(1,*) z(k),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo
 			close(1)
 		endif
 		 		
 		if( z(locnz)==Lz ) then
 			i = 1; j = ny/2 + 1
 			open(1,file='output/debug_data/ustar_of_z_T')
 				do k=1,locnz
 					write(1,*) z(k),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo
 			close(1)
 		endif
 	endif

	
	!--------------------------------------------------------------------------------
	! add psi functions to ustar vector to get u_star_hat, v_star_hat and w_star_hat
	! overwrite ustar, vstar and wstar
	!--------------------------------------------------------------------------------
	if( locnx > 1 .and. .NOT. x_periodic ) ustar(:,:,:) = ustar(:,:,:) + tmpY(:,:,:,4)
	if( ny    > 1 .and. .NOT. y_periodic ) vstar(:,:,:) = vstar(:,:,:) + tmpY(:,:,:,5)
	if( locnz > 1 .and. .NOT. z_periodic ) wstar(:,:,:) = wstar(:,:,:) + tmpY(:,:,:,6)
	
	
	if( istep==istart ) then
		if(myid==0 .and. ny > 1) then
			i=1 ; k=locnz/2
 			open(1,file='output/debug_data/ustar_hat_of_y')
 				do j=1,ny
 					write(1,*) y(j),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo	
 			close(1)
 		endif
 		
 		if(myid==0 .and. locnx > 1) then
			j=1 ; k=locnz/2
 			open(1,file='output/debug_data/ustar_hat_of_x')
 				do i=1,locnx
 					write(1,*) x(i),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo	
 			close(1)
 		endif
 	
 		
 		if( z(1)==0.d0 ) then
 			i = 1; j = ny/2 + 1
 			open(1,file='output/debug_data/ustar_hat_of_z_B')
 				do k=1,locnz
 					write(1,*) z(k),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo
 			close(1)
 		endif
 		 		
 		if( z(locnz)==Lz ) then
 			i = 1; j = ny/2 + 1
 			open(1,file='output/debug_data/ustar_hat_of_z_T')
 				do k=1,locnz
 					write(1,*) z(k),ustar(j,i,k),vstar(j,i,k),wstar(j,i,k)
 				enddo
 			close(1)
 		endif
 	endif

 return
end subroutine update_ustar


subroutine poisson_solver(test_flag)
	!--------------------------------------------------------------- 
	!
	!    grad^2  phi = (1/dt) div ( u*^ )
	!      + homogeneous Neumann conditions
	!
	! After 2d xy transforms we want to solve 1d problems in z
	! 
	!    phi_zz(z,kx,ky) = rhs(z,kx,ky)
	!
	! This routine takes the xy transforms and sets up and calls
	! the 1d z solver. Cosine transforms used in all 3 directions.
	!
	! The rhs =  (1/dt) * div (u*)
	!
	!---------------------------------------------------------------
	use mpi_params,                only: myid
 	use intermediate_variables,    only: phi,div_u,tmpY,tmpZ
 	use independent_variables,     only: dt,nx,ny,nz,y,Lx,Ly,Lz,zg=>z,x_periodic,y_periodic
 	use decomposition_params
 	use differentiation_params,    only: kx,ky,kxfilter,kyfilter
 	use etc,                       only: istep,istart
	implicit none 
	
 	character(len=80),save                      :: exp_type(2)
 	integer,parameter                           :: FOR=1, INV=-1
 	integer                                     :: i,j,k,ig,jg
 	integer, save                               :: locnx,locnz
 	real(kind=8),allocatable,save               :: x(:),z(:)
 	real(kind=8),allocatable                    :: ans(:),rhs(:)
 	real(kind=8)                                :: pi,kval,exact,diff,tol=1.d-10
 	logical                                     :: test_flag
	logical,save                                :: first_entry=.TRUE.
	

	if( first_entry ) then
		!---------------------------------------------------
		!  set expansion type in x & y directions
		!---------------------------------------------------
		exp_type(:) = 'cos'
		if( x_periodic ) exp_type(1)='fourier'
		if( y_periodic ) exp_type(2)='fourier'
		
		locnx = array_size(JDIM,YBLOCK,myid)   ! locnx in YBLOCK
		locnz = array_size(KDIM,YBLOCK,myid)   ! locnz in YBLOCK
		
		!--------------------
		! local x,z vals
		!--------------------
		allocate( x(locnx), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)   ! local x values in YBLOCK
		call get_my_zvals(z,YBLOCK,myid)   ! local z values in YBLOCK
		
		allocate( rhs(nz),ans(nz) )
		!-------------------------------------------------------------------------
		!  TEST THE 1D SOLVER ENGINE
		!------------------------------------------------------------
		!  solve phi'' = rhs  ; phi'=0 at z=0,Lz
		!   phi = cos(2.*pi*z/Lz)  phi_zz = -(2.*pi/Lz)^2 cos(pi*z/Lz)
		!   to do this test, call routine with kx=ky=kval=0.d0
		!-------------------------------------------------------------------------
			pi = 4.d0*atan(1.d0)
			do k=1,nz
				rhs(k) = -(2.*pi/Lz)**2 * cos( 2.*pi*zg(k)/Lz )
			enddo
			kval=0.d0
			call cos_z_solver_phi(ans,rhs,nz,kval,kval)
			do k=1,nz
				exact = cos( 2.*pi*zg(k)/Lz )
				diff = abs(exact-ans(k))
				!write(0,*) ans(k),exact
				if( diff .gt. tol ) stop 'problem in cos_z_solver test problem'
			enddo
			if(myid==0) write(0,*) '...  Analytical test of 1d cos_z_solver successful   ...  tol=1.d-10 '
		deallocate(rhs,ans)
		
		first_entry=.FALSE.
	endif
	
	!-----------------------------------------------------------------------
	!  array/variable map
	!-----------------------------------------------------------------------
	!	rhs => tmpY(:,:,:,1)
	!	rhs_hat => tmpY(:,:,:,2)           ! to store xy transform of rhs
	!	rhs_ZBLOCK  => tmpZ(:,:,:,1)       ! ZBLOCK version of rhs_hat
	!	soln_ZBLOCK => tmpZ(:,:,:,2)       ! ZBLOCK storage format of soln
	!-----------------------------------------------------------------------
 
	!---------------------------------------------------
	!  set the rhs in YBLOCK storage format
	!---------------------------------------------------
	tmpY(:,:,:,1) = (1./dt) * div_u(:,:,:)

	!-------------------------------------------------------------
	! transform the entire 3d rhs array in both x & y directions
	! (transform_xy uses tmpY(:,:,:,6) as wrk space)
	!-------------------------------------------------------------
	call transform_xy(tmpY(:,:,:,1),tmpY(:,:,:,2),FOR,exp_type)

	!--------------------------------------------------------
	! transpose the transformed rhs data to ZBLOCK format,
	! store the result in rhs_ZBLOCK
	!--------------------------------------------------------
	call yblock_2_zblock(tmpY(:,:,:,2),tmpZ(:,:,:,1))
 
	!----------------------------------------------------
	! loop over kx,ky values and do solves
	!----------------------------------------------------
 	do i=1,array_size(JDIM,ZBLOCK,myid)                           ! local x/kx indices in ZBLOCK
 		ig = global_x_indices(START,ZBLOCK,myid) + i - 1          ! global i index
 		do j=1,array_size(KDIM,ZBLOCK,myid)                       ! local y/ky indices in ZBLOCK
 			jg = global_y_indices(START,ZBLOCK,myid) + j - 1      ! global j index
			!----------------------------------------------------
			! do the solve...
			!----------------------------------------------------
			call cos_z_solver_phi(tmpZ(1,i,j,2),tmpZ(1,i,j,1),nz,kx(ig),ky(jg))  ! (soln,rhs,nz,kx,ky)
			tmpZ(:,i,j,2) = tmpZ(:,i,j,2)*kxfilter(ig)*kyfilter(jg)              ! filter in kx,ky
		enddo
	enddo
 
	!------------------------------------------------------------------
	! transpose soln to YBLOCK format, storing result in tmpY(:,:,:,1)
	!------------------------------------------------------------------
	call zblock_2_yblock(tmpZ(1,1,1,2),tmpY(1,1,1,1))
 
	!-------------------------------------------------------------------
	! inverse transform the soln in x and y, store the result in phi 
	!-------------------------------------------------------------------
	call transform_xy(tmpY(1,1,1,1),phi(1,1,1),INV,exp_type)
	
	
	if( istep==istart .and. .not. test_flag) then
		if(myid==0) then
			open(1,file='output/debug_data/phi_of_y')
				k=locnz/2 ; i=1
				do j=1,ny
					write(1,*) y(j),phi(j,i,k),tmpY(j,i,k,1)   ! phi,rhs
				enddo
			close(1)
		endif		
		if( z(1)==0.d0 ) then
			open(1,file='output/debug_data/phi_of_z_B ')
				j=ny/2 + 1 ; i=1
				do k=1,locnz
					write(1,*) z(k),phi(j,i,k),tmpY(j,i,k,1)   ! phi,rhs
				enddo
			close(1)
		endif
		if( z(locnz)==Lz ) then
			open(1,file='output/debug_data/phi_of_z_T ')
				j=ny/2 + 1 ; i=1
				do k=1,locnz
					write(1,*) z(k),phi(j,i,k),tmpY(j,i,k,1)   ! phi,rhs
				enddo
			close(1)
		endif
	endif
 
 return
end subroutine poisson_solver


subroutine test_poisson_solver
!-------------------------------------------------------------------------
!
!     grad^2 phi = rhs   w/ homogeneous Neumann BCs
!
!        to do the test, dt*rhs must be stored in div_u(:,:,:)
!        the solution is stored in phi(:,:,:)
!
!        Let phi = cos(kx)*cos(ly)*cos(mz)
!          phi_xx = -k^2 * phi   
!          phi_yy = -l^2 * phi  
!          phi_zz = -m^2 * phi
!
!          grad^2 phi = -(k^2+l^2+m^2)*cos(kx)*cos(ly)*cos(mz)
!-------------------------------------------------------------------------
	use mpi_params,              only: myid,comm,ierr
	use independent_variables,   only: Lx,Ly,Lz,dt,x_periodic,y_periodic,z_periodic
	use intermediate_variables,  only: tmpY,div_u,phi
	use decomposition_params
	use etc
	implicit none
	real(kind=8), allocatable      :: x(:),y(:),z(:) ! local portions
	integer                        :: i,j,k,locnx,ny,locnz
	real(kind=8)                   :: pi,kk,ll,mm,diff,tol=1.d-10
	logical, save                  :: first_entry=.TRUE., test_flag=.TRUE.
	
	if( first_entry) then
		locnx = array_size(JDIM,YBLOCK,myid)
		ny    = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		
		!--------------------
		! local x,y,z vals
		!--------------------
		allocate( x(locnx), y(ny), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)
		
		pi = 4.d0*atan(1.d0)
		kk = 6.d0 * pi/Lx
		ll = 0.d0 * pi/Ly
		mm = 1.d0 * pi/Lz
		
		! make sure test is periodic if problem is periodic
		if( z_periodic ) mm = 2.d0 * pi/Lz
				
		first_entry = .FALSE.
	endif
	
	!---------------------------------------------------
	! set the rhs: store dt*rhs in div_u
	! (this is poisson solve uses 1/dt * div_u as rhs)
	!---------------------------------------------------
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				div_u(j,i,k) = -(kk**2 + ll**2 + mm**2)*cos(kk*x(i))*cos(ll*y(j))*cos(mm*z(k))
				div_u(j,i,k) = dt*div_u(j,i,k)
				tmpY(j,i,k,3) = cos(kk*x(i))*cos(ll*y(j))*cos(mm*z(k))   ! analytical answer
			enddo
		enddo
	enddo
	
	call poisson_solver(test_flag)    !  the computed solution is stored in phi, routine uses tmpY(:,:,:,1-2)
	
	!--------------------------------
	! check the computed solution
	!--------------------------------
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				diff = abs( phi(j,i,k) - tmpY(j,i,k,3) )
				!if( myid==0 .and. j==1 .and. k==1 ) write(0,*) z(k),phi(j,i,k),tmpY(j,i,k,3),diff 
				if( diff > tol ) stop ' test of poisson_solver failed '
			enddo
		enddo
	enddo
	
	if(myid==0) then
  		message = '...  Analytical test of 3d poisson_solver successful ...  tol=1.d-10'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	
 	div_u = 0.d0
 	phi = 0.d0


 return
end subroutine test_poisson_solver

subroutine cos_z_solver_phi(phi,f,nz,kx,ky)
!------------------------------------------------!
!                                                !
!     solve:    phi_zz(z;kx,ky) = f(z;kx,ky)     !
!                                                !
!               phi_z(0)=phi_z(Lz) = 0           !    
!               z in [0,L]  global length nz     !
!                                                !
!------------------------------------------------!
	use differentiation_params, only: kz,kzfilter,cos_plan,fourier_plan
	use independent_variables,  only: z_periodic
	implicit none
	integer, intent(in)            :: nz
	real(kind=8)                   :: phi(nz),f(nz),kx,ky
	real(kind=8),allocatable,save  :: wrk(:)
 	real(kind=8), save             :: xnorm
 	integer, save                  :: kk, kz_nyq
 	integer(kind=8),save           :: plans(2)
 	integer                        :: k
 	real(kind=8)                   :: xx
 	logical,save                   :: first_entry=.TRUE.
 
	if( first_entry ) then
		allocate( wrk(nz) )
		wrk(:) = 0.d0
		if( z_periodic ) then
			plans(1:2) = fourier_plan(3,1:2)
			xnorm = 1.d0/nz
			kz_nyq = nz/2 + 1             ! location of nyquist wavenumber
			kk=1                          ! fourier never trims endpts
		else
			plans(1:2) = cos_plan(3)
			xnorm = 1.d0/(2.d0*(nz-1.d0))
			kz_nyq = nz                   ! location of nyquist wavenumber
			kk=1                          ! don't trim endpts, cos expansion
		endif
		
				
		first_entry=.FALSE.
	endif
	
	phi(1:nz) = 0.d0

	!---------------------------------------
	!  take forward transform of the rhs f
	!---------------------------------------
	call dfftw_execute_r2r(plans(1),f(kk),wrk(kk)) 

	!---------------------------------------
	!  divide by -(kx^2+ky^2+kz^2)
	!---------------------------------------
	do k=2,nz
		xx = -( kx**2 + ky**2 + kz(k)**2 )
		wrk(k) = xnorm*wrk(k)/xx
	enddo
 
	!-----------------------------------------
	! zero out the kx=ky=kz=0 mean component
	! this is the compatability condition for
	! homogeneous Neumann: mean of rhs=0,
	! this just avoids the divide by zero
	!-----------------------------------------
	if( kx==0.d0 .and. ky==0.d0 ) then
		wrk(1)=0.d0
	else                             ! kz=0 but kx^2 + ky^2 .NE. 0
		xx = -( kx**2 + ky**2 )
		wrk(1) = xnorm*wrk(1)/xx
	endif
 
	!--------------------------------------------
	! Apply kzfilter
 	!--------------------------------------------
 	wrk(:) = wrk(:)*kzfilter(:)
 	
 	!--------------------------------------------
	! Nyquist always zero
 	!--------------------------------------------
 	wrk(kz_nyq) = 0.d0
 
	!---------------------------------------
	!  take inverse transform of result
	!---------------------------------------
	call dfftw_execute_r2r(plans(2),wrk(kk),phi(kk))
 
 return
end subroutine cos_z_solver_phi





