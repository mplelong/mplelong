subroutine diffusion_coeffs
	!-----------------------------------------------------------------------------------------
	!   Given user input, set up the diffusion coefficients. Regardless of the
	!   type of diffusion (Laplacian or high_order_operators) the diffusion
	!   operators take the form (for example)
	!
	!   -( nu_star(1)*kx^(2p(1)) + nu_star(2)*kx^(2p(2)) + nu_star(3)*kx^(2p(3)) ) u^(kx,ky,kz)
	!
	!   when high_order_operators is false, the nu_star values are all set to nu and all the p values to 1
	!   similar treatment for kappa_star for scalars 1 and (possibly) 2
	!
	!   when  high_order_operators is true, user specified diffusion times (at grid scales) and p values
	!   are combined with the nyquist wavenumbers to compute the nu_star and kappa_star values
	!-----------------------------------------------------------------------------------------
	use mpi_params,              only: myid,comm,ierr
	use dimensional_scales,      only: nu,kappa, nu_star,kappa_star,T_diff  ! nu_star(3)m kappa_star(3,2) x,y,z
	use dimensionless_params,    only: p                                    ! integer (3) x,y,z, half orders
	use differentiation_params,  only: kx,ky,kz                             ! wavenumbers already computed
	use independent_variables,   only: nx,ny,nz                             ! to access the largest wavenumber values
	use methods_params,          only: high_order_operators,do_second_scalar
	implicit none
	real(kind=8)                    :: denom
	integer                         :: idir
	
	
	if( .not. high_order_operators ) then
	
		p(:) = 1                          ! 1/2 order for Laplacian diffusion
		nu_star(:) = nu                   ! [m2/s] 
		kappa_star(:,1) = kappa(1)        ! [m2/s]
		if(do_second_scalar) then
			kappa_star(:,2) = kappa(2)    ! [m2/s]
		endif
		
		! these time scales not actually used, just calculated to be able to write the values during initialization
		T_diff(:) = 1.d23   ! ~ infinite
		if( nx > 1 ) then
			T_diff(1) = 1.d0/(kx(nx)**(2*p(1))*nu_star(1))    ! calculate time scale for momentum diffusion at x nyquist scale
		endif
		if( ny > 1 ) then
			T_diff(2) = 1.d0/(ky(ny)**(2*p(2))*nu_star(2))    ! calculate time scale for momentum diffusion at x nyquist scale
		endif
		if( nz > 1 ) then
			T_diff(3) = 1.d0/(kz(nz)**(2*p(3))*nu_star(3))    ! calculate time scale for momentum diffusion at x nyquist scale
		endif
		
	elseif( high_order_operators) then
	
		nu_star(:) = 0.d0
		kappa_star(:,:) = 0.d0
		
		! choose nu* such that T_diff = 1./(k_max^(2p) * nu*)
		idir = 1
		if( nx > 1 ) then
			denom = kx(nx)**(2*p(idir)) * T_diff(idir)
			nu_star(idir) = 1.d0/denom
		endif
	
		idir = 2
		if( ny > 1 ) then
			denom = ky(ny)**(2*p(idir)) * T_diff(idir)
			nu_star(idir) = 1.d0/denom
		endif
	
		idir = 3
		if( nz > 1 ) then
			denom = kz(nz)**(2*p(idir)) * T_diff(idir)
			nu_star(idir) = 1.d0/denom
		endif
	
		!  set high order scalar and momentum diffusion coefficients equal to one another
		kappa_star(:,1) = nu_star(:)
		
		!  set both scalar diffusion coefficients equal to one another
		if(do_second_scalar) then
			kappa_star(:,2) = kappa_star(:,1)
		endif
		
	endif	
	
	if( myid==0 ) then
		open(1,file='output/debug_data/diffusion_coeffs')
			idir=1
			write(1,*) kx(nx), T_diff(idir), p(idir), nu_star(idir), kappa_star(idir,1)
			idir=2
			write(1,*) ky(ny), T_diff(idir), p(idir), nu_star(idir), kappa_star(idir,1)
			idir=3
			write(1,*) kz(nz), T_diff(idir), p(idir), nu_star(idir), kappa_star(idir,1)
		close(1)
	endif
	

	call mpi_barrier(comm,ierr)	
 return
end subroutine diffusion_coeffs	
	
subroutine diffuse
	!-----------------------------------------------------------------------------------
	! Routine to integrate the diffusive terms analytically using an integrating
	! factor. After the explicit time integration and the pressure projection step,
	! we can write
	!               u_np1 = u** + int diffussive term dt
	!                   where u** is the result after explicit time integration
	!                             and pressure projection
	!
	! Transforming to wavenumber space, this is equivalent to integrating
	!              d/dt u_hat = -nu_* k^2p u_hat
	! for one time step with the initial condition u_hat = u**_hat
	!
	! The exact solution to this eqn. is: u_hat = u**_hat * exp( - nu_* k^2p dt )
	! and u_np1 is then just the inverse transform of u_hat.
	!
	! This routine implements this logic in 3 dimensions for velocity and scalars.
	!
	!-----------------------------------------------------------------------------------
	use mpi_params,                   only: myid
 	use decomposition_params
 	use independent_variables,        only: dt,ny,nz              ! need global nz to pass to z_diffusion
 	use dependent_variables,          only: u,v,v,w,s1,s2,s1_bar,s2_bar
	use intermediate_variables,       only: tmpX,tmpY,tmpZ
	use differentiation_params,       only: kx,ky,kxfilter,kyfilter
	use methods_params,               only: do_second_scalar,do_sponging
	use dimensionless_params,         only: p                     ! 1/2 order of diffusion operators, x,y,z
	use dimensional_scales,           only: nu_star, kappa_star   ! nu/kappa_star(3) x,y and z dirs
 	implicit none
 	integer, save                        :: nvars=4,locnx,locnz
 	integer                              :: i,j,k,ig,id,fid,dir
 	real(kind=8)                         :: xx,yy
 	real(kind=8),allocatable,save        :: diff_factor(:,:,:)
 	character(len=80),save               :: exp_type='cos'
 	real(kind=8), external               :: myexp
 	logical, save                        :: first_entry=.TRUE.
	
	if( first_entry) then
	
		if( do_second_scalar )  nvars = 5
		
		locnx = array_size(JDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		allocate( diff_factor(ny,locnx,3) )
		
		do i=1,locnx
			ig = global_x_indices(START,YBLOCK,myid) + i - 1
			do j=1,ny
			
				xx = nu_star(1)*kx(ig)**(2*p(1))
				yy = nu_star(2)*ky(j)**(2*p(2))
				diff_factor(j,i,1) = myexp(-(xx+yy)*dt)      ! for u,v,w
				
				id=1  ! scalar 1
				xx = kappa_star(1,id)*kx(ig)**(2*p(1))
				yy = kappa_star(2,id)*ky(j)**(2*p(2))
				diff_factor(j,i,2) = myexp(-(xx+yy)*dt)      ! for s1
				
				if( do_second_scalar ) then
					id=2  ! scalar 2
					xx = kappa_star(1,id)*kx(ig)**(2*p(1))
					yy = kappa_star(2,id)*ky(j)**(2*p(2))
					diff_factor(j,i,3) = myexp(-(xx+yy)*dt)  ! for s2
				endif
				
				! apply wavenumber filtering to high wavenumbers
				diff_factor(j,i,:) = diff_factor(j,i,:)*kxfilter(ig)*kyfilter(j)
								
			enddo
		enddo
				
		first_entry = .FALSE.
	endif
	
	do id=1,nvars
		!---------------------------------------------------------------------------------
		!  transform variable "id" in x and y directions, result is is returned as YBLOCK
		!---------------------------------------------------------------------------------
		dir = 1   ! FOR transforms (x,y) -> (kx,ky)
		if(id==1) then
			call transform_xy(u,tmpY,dir,exp_type)
			fid = 1
		elseif(id==2) then 
			call transform_xy(v,tmpY,dir,exp_type)
			fid = 1
		elseif(id==3) then 
			call transform_xy(w,tmpY,dir,exp_type)
			fid = 1
		elseif(id==4) then 
			call transform_xy(s1,tmpY,dir,exp_type)
			fid = 2
		elseif(id==5) then
			call transform_xy(s2,tmpY,dir,exp_type)
			fid = 3
		endif
	
		!------------------------------------------------------------------------
		!  multiply by exp(-nu_* (kx^2p + ky^2p)*dt), correct version from fid
		!------------------------------------------------------------------------
		do k=1,locnz
			do i=1,locnx
				do j=1,ny
					tmpY(j,i,k,1) = tmpY(j,i,k,1)*diff_factor(j,i,fid)
				enddo
			enddo
		enddo
		
	
		!------------------------------------------------------------------------
		!  transpose to ZBLOCK(z,x,y)
		!------------------------------------------------------------------------
		call YBLOCK_2_ZBLOCK(tmpY,tmpZ)
	
		do i=1,array_size(JDIM,ZBLOCK,myid)
			do j=1,array_size(KDIM,ZBLOCK,myid)
				!-------------------------------------------------------------------
				! at each kx,ky, transform in z,  apply z factor, inverse transform
				!-------------------------------------------------------------------
				call z_diffusion(tmpZ(1,i,j,1),tmpZ(1,i,j,2),nz,id)
			enddo
		enddo
	
		!------------------------------------------------------------------------
		!  transpose back to YBLOCK, storing result in tmpY(1)
		!------------------------------------------------------------------------
		call ZBLOCK_2_YBLOCK(tmpZ(1,1,1,2),tmpY(1,1,1,1))
	
	
		!------------------------------------------------------------------------
		!  inverse transform in x and y, overwrite dependent variables with
		!  their diffusively damped updates. 
		!------------------------------------------------------------------------
		dir = -1  ! INV transforms (kx,ky) -> (x,y)
		if(id==1) call transform_xy(tmpY,u,dir,exp_type)
		if(id==2) call transform_xy(tmpY,v,dir,exp_type)
		if(id==3) call transform_xy(tmpY,w,dir,exp_type)
		if(id==4) call transform_xy(tmpY,s1,dir,exp_type)
		if(id==5) call transform_xy(tmpY,s2,dir,exp_type)
		
				
	enddo   ! end id loop over nvars
 return		
end subroutine diffuse




subroutine z_diffusion(f,ans,nz,id)
!--------------------------------------------------!
!                                                  !
!     Given f(z) take z transform,  apply the      !
!     wavenumber space damping factor, and         !
!     inverse transform, returning the damped      !
!     field in the array ans. The correct damping  !
!     factor is selected using the variable id.    !
!                                                  !
!     This method imposes adiabatic boundary       !
!     conditions  ==> if diffusive fluxes are      !
!     needed then additional source terms must     !
!     instead                                      !
!--------------------------------------------------!
	use differentiation_params, only: kz,kzfilter,cos_plan  ! cos_plan(3) in x,y,z dirs
	use dimensionless_params,   only: p                     ! 1/2 order of diffusion operators, x,y,z
	use dimensional_scales,     only: nu_star, kappa_star   ! nu_star(3)/kappa_star(3,2) x,y and z dirs, scalar 1/2
	use independent_variables,  only: dt
	implicit none
	integer, intent(in)            :: nz,id
	real(kind=8)                   :: f(nz),ans(nz)
	real(kind=8),allocatable,save  :: wrk(:),diff_factor(:,:)
 	real(kind=8),save              :: xnorm
 	integer,save                   :: kk=1  ! transform full data array, i.e. start at position 1
 	integer(kind=8),save           :: plans(2)
 	integer                        :: k,fid
 	real(kind=8)                   :: zz
 	real(kind=8), external         :: myexp
 	logical,save                   :: first_entry=.TRUE.
 
	if( first_entry ) then
		allocate( wrk(nz), diff_factor(nz,3) )
		wrk(:) = 0.d0
		diff_factor(:,:) = 0.d0
		
		plans(1:2) = cos_plan(3)
		xnorm = 1.d0/(2.d0*(nz-1.d0))
		
		!-----------------------------------------------------------
		! save the damping factors for the z direction (index 3)
		!-----------------------------------------------------------
		do k=1,nz
			zz = nu_star(3)*kz(k)**(2*p(3))
			diff_factor(k,1) = myexp(-zz*dt)        ! for u,v,w
			
			zz = kappa_star(3,1)*kz(k)**(2*p(3))
			diff_factor(k,2) = myexp(-zz*dt)        ! for s1
			
			zz = kappa_star(3,2)*kz(k)**(2*p(3))
			diff_factor(k,3) = myexp(-zz*dt)        ! for s2
			
			! apply wavenumber filtering to high wavenumbers
			diff_factor(k,:) = diff_factor(k,:)*kzfilter(k)
		enddo
				
		first_entry=.FALSE.
	endif
	
	ans(1:nz) = 0.d0      ! why should this be important???
	
	if( id <= 3 ) then
		fid = 1    ! use diff_factor for u,v,w (nu_star(:))
	elseif( id==4 ) then
		fid = 2    ! use diff_factor for s1 (kappa_star(:,1))
	elseif( id==5 ) then
		fid = 3    ! use diff_factor for s1 (kappa_star(:,2))
	endif

	!----------------------------------------------
	!  take forward transform of f(z)
	!        wrk <-- f_hat(kz)
	!----------------------------------------------
	call dfftw_execute_r2r(plans(1),f(kk),wrk(kk)) 

	!----------------------------------------------
	!  multiply by the damping factor and
	!  normalize the result, overwriting wrk(:)
	!----------------------------------------------
	do k=1,nz
		wrk(k) = xnorm * diff_factor(k,fid) * wrk(k)
	enddo
  	
 	!--------------------------------------------
	! Nyquist always explicitly set to zero
 	!--------------------------------------------
 	wrk(nz) = 0.d0
 
	!---------------------------------------------------------
	!  take inverse transform of the damped, normalized data,
	!  storing the result in ans(:)
	!---------------------------------------------------------
	call dfftw_execute_r2r(plans(2),wrk(kk),ans(kk))
 
 return
end subroutine z_diffusion

subroutine test_z_diffusion
	use mpi_params,                   only: myid
	use independent_variables,        only: z,nz,Lz,dt
	use dimensionless_params,         only: p           ! 1/2 order of diffusion operators, x,y,z
	use dimensional_scales,           only: nu_star     ! nu_star(3) x,y and z dirs
	use etc
	implicit none
	real(kind=8), allocatable            :: phi(:),phi_exact(:),phi_computed(:)
	real(kind=8)                         :: t0,tf,kappa,z0,dz,original_nu_star,arg
	real(kind=8)                         :: error, tol = 1.d-12
	integer                              :: k,original_p,id
	real(kind=8), external               :: myerf
	
	allocate( phi(nz),phi_exact(nz),phi_computed(nz) )
	phi = 0.d0
	phi_exact = 0.d0
	phi_computed = 0.d0
	z0 = Lz/2.d0
	t0 = 10*dt                      ! initial field is smooth
	tf = t0 + dt                    ! integrate exactly 1 time step (analytically)
	
	original_nu_star = nu_star(3)   ! z dim value
	original_p = p(3)               ! z dim value
	
	dz = z(2)-z(1)
	nu_star(3) = 0.5d0 * (dz**2)/dt   ! choose a reasonable test value for the diffusion coeff
	p(3) = 1                          ! 1/2 order; test Laplacian diffusion
	
	!---------------------------------------------------
	! set the initial condition:  exact soln at t=t0
	! and the exact solution  at tf, 1 time step later
	!---------------------------------------------------
	do k=1,nz
		arg = (z(k)-z0) / sqrt(4.d0*nu_star(3)*t0)
		phi(k) = 0.5d0*( 1.d0 + myerf(arg) )
		
		arg = (z(k)-z0) / sqrt(4.d0*nu_star(3)*tf)
		phi_exact(k) = 0.5d0*( 1.d0 + myerf(arg) )		
	enddo
	
	! integrate forward 1 dt time step
	id = 1        ! i.e. have z_diffusion use nu_star and not kappa_star
	call z_diffusion(phi,phi_computed,nz,id)
	
	! check the solution
	do k=1,nz
		error = abs( phi_computed(k) - phi_exact(k) )
		!write(0,*) phi_computed(k),phi_exact(k),error
		if( myid==0 .and. error > tol ) stop ' test of z_diffusion failed, tol=1.d-12 '
	enddo
	
	!  replace original values
	nu_star(3) = original_nu_star
	p(3) = original_p
	
	! deallocate tmp arrays
	deallocate( phi, phi_exact, phi_computed )
	
	if(myid==0) then
  		message = '...  Analytical test of 1d z_diffusion successful    ...  tol=1.d-12'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 return	
end subroutine test_z_diffusion


!to compute grad2 u  as div( grad u )
!	(1) grad u    		2 transposes ( Y <-> X ) ( Y <-> Z )   ( 1 transpose = FOR and INV )
!						3 transforms (x,y,z)                   ( 1 transform = FOR and INV )
!						3 arrays to hold u_x, u_y and u_z
!					
!	(2) div( grad u)	2 transposes ( Y <-> X ) ( Y <-> Z )
!						3 transforms (x,y,z)
!						
!	total:   4 transposes 6 transforms 3 arrays


	
!to apply diffusion to u in wavenumber space, i.e. use integrating factor:

!	(1) xy_transform:   y transform, Y->X, x transform, X->Y:  	1 transpose ( FOR and INV)
!																2  1/2 transforms (FOR only)
!																
!	(2) transpose Y->Z, z transform, apply factor, inverse transform, Z->Y
!																1 transpose (FOR and INV)
!																1 transform (FOR and INV)
																
!    (3) inv xy transform: y transform, Y->X, x transform, X->Y: 1 transpose ( FOR and INV)
!																2  1/2 transforms (FOR only)
																
!	total:     3 transposes, 3 transforms      THIS IS MUCH LESS WORK



