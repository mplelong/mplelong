SUBROUTINE deriv_BC(f,df,n,dir,method,debug)
!--------------------------------------------------------------------------------------------
!
!    differentiate a 1d data array of length n, no symmetries assumed
!
!    compute deriv f wrt 'dir' given f at equally spaced grid points in [0,L] using the
!    Bernoulli-cosine approach to deal with the discontinuous derivatives at 0,L
!    for the even extension of f into [0,2L)
!
!    dir = 1,2,3  for x, y and z directions respectively
!    tmp has to be big enough to contain an (n,4) array
!
!    frac is the fraction of the wavenumber range to filter in cos transform differentiation
!         frac<0 ==> 2/3 rule, frac0 ==> no filtering
!    Q is the highest odd order of the U_n terms in the singular series expansions
!         (Q+1)/2 is the number of terms in the expansions
!    LU = [lu_A,lu_B,piv_A,piv_B] the prefactored matrices constraining the expansion coeffs
!
!    if Q=0, just do standard (filtered) term x term cosine differentiation
!
!    data needed to carry out the differentiation obtained from modules
!---------------------------------------------------------------------------------------------
	use mpi_params,                          only: myid,comm,ierr
	use independent_variables,               only: x,y,z,nx,ny,nz,Lx,Ly,Lz
	use differentiation_params,              only: LU_x,ipiv_x,LU_y,ipiv_y,LU_z,ipiv_z, &
	                                               U_x,dU_x,U_y,dU_y,U_z,dU_z,Q
	use fourier_differentiation_tools,       only: differentiate_fcs
	implicit none
	logical, intent(in)                         :: debug
	real(kind=8), intent(in)                    :: f(n)
	real(kind=8), intent(inout)                 :: df(n)
	integer, intent(in)                         :: n,dir
	character(len=80)                           :: method,BCmethod
	integer, parameter                          :: order=1
	real(kind=8)                                :: x0,L  
	integer                                     :: j,K,nmax
	integer, save                               :: M
	real(kind=8), allocatable, save             :: LU_0(:,:), LU_L(:,:)
	real(kind=8), allocatable, save             :: rhs(:),A(:),B(:),tmp(:,:)
	integer, allocatable, save                  :: ipiv_0(:), ipiv_L(:)
	real(kind=8), allocatable, save             :: U(:,:,:), dU(:,:,:)
	logical, save                               :: first_entry=.TRUE.
		
	if(first_entry) then
		M = (Q+1)/2
		allocate( LU_0(M,M), LU_L(M,M), ipiv_0(M), ipiv_L(M), rhs(M), A(M), B(M) )
		LU_0=0.d0 ; LU_L=0.d0; rhs=0.d0; A=0.d0; B=0.d0
		
		nmax = MAXVAL( (/nx,ny,nz/) )
		allocate( tmp(nmax,4) )
		tmp = 0.d0
		
		allocate( U(nmax,M,2), dU(nmax,M,2) )  ! big enough for all 3 directions		
		first_entry=.FALSE.
	endif
	
		
	if( n==1 ) then
		df=0.d0
		return
	endif
	
	!-----------------------------------------------------------
	!  straight Fourier method, 'fourier', 'cos' or 'sin'
	!-----------------------------------------------------------	
	if( method=='fourier' .or. method=='cos' .or. method=='sin' ) then
		call differentiate_fcs(f,df,n,dir,method,order)
		return	
	endif
	
	if( dir==3 .and. debug .and. myid==0 ) then
		write(0,*) '                                    debug ON for deriv_BC dir, n: ',dir,n
	endif
	
	!-----------------------------------------------------------
	!  for Bernoulli-Cosine, get stuff for correct dir
	!-----------------------------------------------------------
	tmp(1:n,:) = 0.d0
	if(dir==1) then
		L=Lx
		LU_0(1:M,1:M)   =   LU_x(1:M,1:M,1)
		LU_L(1:M,1:M)   =   LU_x(1:M,1:M,2)
		ipiv_0(1:M)     = ipiv_x(1:M,1)
		ipiv_L(1:M)     = ipiv_x(1:M,2)
		U(1:nx,1:M,1:2)   =  U_x(1:nx,1:M,1:2)    ! U_x(x,j,K) is the jth basis function in x direction expanded about 0 or Lx (K)
		dU(1:nx,1:M,1:2)  = dU_x(1:nx,1:M,1:2)    ! dU_x(x,j,K) is the deriv of jth basis function in x direction expanded about 0 or Lx (K)
	elseif(dir==2) then
		L=Ly
		LU_0(1:M,1:M)   =   LU_y(1:M,1:M,1)
		LU_L(1:M,1:M)   =   LU_y(1:M,1:M,2)
		ipiv_0(1:M)     = ipiv_y(1:M,1)
		ipiv_L(1:M)     = ipiv_y(1:M,2)
		U(1:ny,1:M,1:2)   =  U_y(1:ny,1:M,1:2)    ! U_y(x,j,K) is the jth basis function in x direction expanded about 0 or Ly (K)
		dU(1:ny,1:M,1:2)  = dU_y(1:ny,1:M,1:2)    ! dU_y(x,j,K) is the deriv of jth basis function in x direction expanded about 0 or Ly (K)
	elseif(dir==3) then
		L=Lz
		LU_0(1:M,1:M)   =   LU_z(1:M,1:M,1)
		LU_L(1:M,1:M)   =   LU_z(1:M,1:M,2)
		ipiv_0(1:M)     = ipiv_z(1:M,1)
		ipiv_L(1:M)     = ipiv_z(1:M,2)
		U(1:nz,1:M,1:2)   =  U_z(1:nz,1:M,1:2)    ! U_z(x,j,K) is the jth basis function in x direction expanded about 0 or Lz (K)
		dU(1:nz,1:M,1:2)  = dU_z(1:nz,1:M,1:2)    ! dU_z(x,j,K) is the deriv of jth basis function in x direction expanded about 0 or Lz (K)
	endif
	
	
			
	!---------------------------------------------------------------
	!  solve the matrix problems for coeffs for 2 series expansions
	!---------------------------------------------------------------
	x0 = 0.d0
	call solve_for_B_expansion_coeffs(LU_0,ipiv_0,f,n,M,x0,rhs,A)
		
	x0 = L
	call solve_for_B_expansion_coeffs(LU_L,ipiv_L,f,n,M,x0,rhs,B)
	
	if( dir==3 .and. debug .and. myid==0 ) then
		open(1,file='output/debug_data/expansion_coeffs')
		do j=1,M
			write(1,*) j,A(j),B(j)
		enddo
		close(1)
	endif
	
	
	!---------------------------------------------------------------
	!  construct the 2 series for x in [0,L], add them to get f_s(x)
	!  also differentiate this series
	!---------------------------------------------------------------
	do j=1,M
		K = 1
		tmp(1:n,1)   = tmp(1:n,1)   + A(j)* U(1:n,j,K)  ! s0     i.e. series expanded about 0
		tmp(1:n,2)   = tmp(1:n,2)   + A(j)*dU(1:n,j,K)  ! s0_x   deriv of series expanded about 0
		
		K = 2	
		tmp(1:n,3)   = tmp(1:n,3)   + B(j)* U(1:n,j,K)  ! sL     i.e. series expanded about L
		tmp(1:n,4)   = tmp(1:n,4)   + B(j)*dU(1:n,j,K)  ! sL_x   deriv of series expanded about L
	enddo
	if( dir==3 .and. debug .and. myid==0 ) then
		open(1,file='output/debug_data/separate_series')
		do j=1,n
			write(1,*) z(j),tmp(j,1),tmp(j,2),tmp(j,3),tmp(j,4)
		enddo
		close(1)
	endif
	
	
	 	
	! combine the A/0 and B/L series;  store in tmp arrays
	tmp(1:n,1) = tmp(1:n,1) + tmp(1:n,3)    ! f_s   = s0   + sL
	tmp(1:n,2) = tmp(1:n,2) + tmp(1:n,4)    ! f_s_x = s0_x + sL_x
	if( dir==3 .and. debug .and. myid==0 ) then
		open(1,file='output/debug_data/combined_series')
		do j=1,n
			write(1,*) z(j),tmp(j,1),tmp(j,2),f(j)
		enddo
		close(1)
	endif

	! extract the smooth, even extendable part that is well approximated by a cosine series
	tmp(1:n,1) = f(1:n) - tmp(1:n,1)  ! i.e. f_Q = f - f_s  should be Q times differentiable when even expanded 
	if( dir==3 .and. debug .and. myid==0 ) then
		open(1,file='output/debug_data/smooth_part')
		do j=1,n
			write(1,*) z(j),tmp(j,1)
		enddo
		close(1)
	endif
	
		
	!---------------------------------------------------------------
	!  use standard cosine transform to differentiate f_Q
	!  store result in tmp(:,3) use tmp(:,4) as work space
	!---------------------------------------------------------------
	BCmethod='cos'
	call differentiate_fcs(tmp(1,1),tmp(1,3),n,dir,method,order)
	
	
		
	df(1:n) = tmp(1:n,3) + tmp(1:n,2)  ! i.e. f_Q_x + f_s_x
	if( dir==3 .and. debug .and. myid==0 ) then
		open(1,file='output/debug_data/BCderivs')
		do j=1,n
			write(1,*) z(j),tmp(j,3),tmp(j,2),df(j)
		enddo
		close(1)
	endif

 return
END SUBROUTINE deriv_BC







!-------------------------------------------------------------------------
SUBROUTINE evaluate_basis_functions
!-------------------------------------------------------------------------
!    Construct the U_n basis functions for n=1,3,5... and their derivatives
!    at each x point for the (Q+1)/2 term series
!    The S_0 series are expended about x=0 while 
!    the S_L series are expanded about x=x[-1]=L
!    Do this for each of the x,y and z directions. The results become available
!    in the module differentiation_params as  
!       U_x(x,n,K) dU_x(x,n,K)   U_y(y,n,K) dU_y(y,n,K)    U_z(z,n,K) dU_z(z,n,K)
!-----------------------------------------------------------------------------------
	use mpi_params,                only: myid,comm,ierr
	use differentiation_params
	use independent_variables,     only: x,y,z,nx,ny,nz,Lx,Ly,Lz
	implicit none
	real(kind=8), external            :: U_n, ddx_U_n
	integer                           :: M,i,j,n
	real(kind=8)                      :: a,b,P
	
	M = (Q+1)/2     ! number of terms in the series expansions
	allocate( U_x(nx,M,2), dU_x(nx,M,2) )
	allocate( U_y(ny,M,2), dU_y(ny,M,2) )
	allocate( U_z(nz,M,2), dU_z(nz,M,2) )
	U_x=0.d0 ; dU_x = 0.d0
	U_y=0.d0 ; dU_y = 0.d0
	U_z=0.d0 ; dU_z = 0.d0
	
	!  x direction
	if( nx > 1 ) then
		a = 0.d0        ! start point
		b = Lx          ! end point 
		P = 2.d0*Lx     ! periodicity length for even extension
		do i=1,nx
			do j=1,M
				n = 2*j-1
				U_x(i,j,1) = U_n(x(i),a,P,n)
				U_x(i,j,2) = U_n(x(i),b,P,n)
			
				dU_x(i,j,1) = ddx_U_n(x(i),a,P,n)
				dU_x(i,j,2) = ddx_U_n(x(i),b,P,n)
			enddo
		enddo
	endif
	
	!  y direction
	if( ny > 1 ) then
		a = 0.d0        ! start point
		b = Ly          ! end point 
		P = 2.d0*Ly     ! periodicity length for even extension
		do i=1,ny
			do j=1,M
				n = 2*j-1
				U_y(i,j,1) = U_n(y(i),a,P,n)
				U_y(i,j,2) = U_n(y(i),b,P,n)
			
				dU_y(i,j,1) = ddx_U_n(y(i),a,P,n)
				dU_y(i,j,2) = ddx_U_n(y(i),b,P,n)
			enddo
		enddo
	endif
	
	!  z direction
	if( nz > 1 ) then
		a = 0.d0        ! start point
		b = Lz          ! end point 
		P = 2.d0*Lz     ! periodicity length for even extension
		do i=1,nz
			do j=1,M
				n = 2*j-1
				U_z(i,j,1) = U_n(z(i),a,P,n)
				U_z(i,j,2) = U_n(z(i),b,P,n)
			
				dU_z(i,j,1) = ddx_U_n(z(i),a,P,n)
				dU_z(i,j,2) = ddx_U_n(z(i),b,P,n)
			enddo
		enddo
	endif
	
	! x
	if(myid==0) then
		n=1
		open(1,file='output/debug_data/x_bases_1')
		do i=1,nx
			write(1,*) x(i),U_x(i,n,1),U_x(i,n,2),dU_x(i,n,1),dU_x(i,n,2)
		enddo
		close(1)
		
		n=2
		open(1,file='output/debug_data/x_bases_2')
		do i=1,nx
			write(1,*) x(i),U_x(i,n,1),U_x(i,n,2),dU_x(i,n,1),dU_x(i,n,2)
		enddo
		close(1)
		
		n=3
		open(1,file='output/debug_data/x_bases_3')
		do i=1,nx
			write(1,*) x(i),U_x(i,n,1),U_x(i,n,2),dU_x(i,n,1),dU_x(i,n,2)
		enddo
		close(1)
		
		if( M>=4 ) then
			n=4
			open(1,file='output/debug_data/x_bases_4')
			do i=1,nx
				write(1,*) x(i),U_x(i,n,1),U_x(i,n,2),dU_x(i,n,1),dU_x(i,n,2)
			enddo
			close(1)
		endif
		
		if( M>=5 ) then
			n=5
			open(1,file='output/debug_data/x_bases_5')
			do i=1,nx
				write(1,*) x(i),U_x(i,n,1),U_x(i,n,2),dU_x(i,n,1),dU_x(i,n,2)
			enddo
			close(1)
		endif
	endif
	
	! y
	if(myid==0) then
		n=1
		open(1,file='output/debug_data/y_bases_1')
		do i=1,ny
			write(1,*) y(i),U_y(i,n,1),U_y(i,n,2),dU_y(i,n,1),dU_y(i,n,2)
		enddo
		close(1)
		
		n=2
		open(1,file='output/debug_data/y_bases_2')
		do i=1,ny
			write(1,*) y(i),U_y(i,n,1),U_y(i,n,2),dU_y(i,n,1),dU_y(i,n,2)
		enddo
		close(1)
		
		n=3
		open(1,file='output/debug_data/y_bases_3')
		do i=1,ny
			write(1,*) y(i),U_y(i,n,1),U_y(i,n,2),dU_y(i,n,1),dU_y(i,n,2)
		enddo
		close(1)
		
		if( M>=4 ) then
			n=4
			open(1,file='output/debug_data/y_bases_4')
			do i=1,ny
				write(1,*) y(i),U_y(i,n,1),U_y(i,n,2),dU_y(i,n,1),dU_y(i,n,2)
			enddo
			close(1)
		endif
		
		if( M>=5 ) then
			n=5
			open(1,file='output/debug_data/y_bases_5')
			do i=1,ny
				write(1,*) y(i),U_y(i,n,1),U_y(i,n,2),dU_y(i,n,1),dU_y(i,n,2)
			enddo
			close(1)
		endif
	endif
	
	! z
	if(myid==0) then
		n=1
		open(1,file='output/debug_data/z_bases_1')
		do i=1,nz
			write(1,*) z(i),U_z(i,n,1),U_z(i,n,2),dU_z(i,n,1),dU_z(i,n,2)
		enddo
		close(1)
		
		n=2
		open(1,file='output/debug_data/z_bases_2')
		do i=1,nz
			write(1,*) z(i),U_z(i,n,1),U_z(i,n,2),dU_z(i,n,1),dU_z(i,n,2)
		enddo
		close(1)
		
		n=3
		open(1,file='output/debug_data/z_bases_3')
		do i=1,nz
			write(1,*) z(i),U_z(i,n,1),U_z(i,n,2),dU_z(i,n,1),dU_z(i,n,2)
		enddo
		close(1)
		
		if( M>=4 ) then
			n=4
			open(1,file='output/debug_data/z_bases_4')
			do i=1,nz
				write(1,*) z(i),U_z(i,n,1),U_z(i,n,2),dU_z(i,n,1),dU_z(i,n,2)
			enddo
			close(1)
		endif
		
		if( M>=5 ) then
			n=5
			open(1,file='output/debug_data/z_bases_5')
			do i=1,nz
				write(1,*) z(i),U_z(i,n,1),U_z(i,n,2),dU_z(i,n,1),dU_z(i,n,2)
			enddo
			close(1)
		endif
	endif
	
 call mpi_barrier(comm,ierr)
 return
END SUBROUTINE evaluate_basis_functions	





!-------------------------------------------------------------------------
SUBROUTINE setup_factor_B_matrices
!-------------------------------------------------------------------------
!    Set up the matrix problem that determines the weights of the U_n(x) functions
!    x is a vector of equally spaced points [0,L] where f is defined.
!    The even extension of f is periodic over length 2L. 
!    x0 is the endpoint around which the U_n functions are expanded 
!    (i.e. x=0 and x=L) for the A and B series expansions. The matrix eqns
!    constrain the Q term expansions to match f at  (Q+1)/2 points near 0 and L.
!
!    The factored coefficient matrices are made available via the
!    module differentiation_params as LU_x(:,:,K), LU_y(:,:,K), LU_z(:,:,K) 
!    and ipiv_x(:,K), ipiv_x(:,K), ipiv_x(:,K)  where K=1 for expansions
!    around the start point and K=1 for expansions around the end point 
!--------------------------------------------------------------------------------------
	use mpi_params,                       only: myid,comm,ierr
	use differentiation_params,           only: Q, LU_x, LU_y, LU_z, ipiv_x, ipiv_y, ipiv_z
	use independent_variables,            only: x,y,z,nx,ny,nz,Lx,Ly,Lz
	use mpi_params,                       only: myid
	implicit none
	integer                                  :: npts,info,K,istart,i,j,n
	real(kind=8)                             :: a,L,P,xval,yval,zval
	real(kind=8), external                   :: U_n
	
	npts = (Q+1)/2                        ! number of terms in expansion, values to match near each endpoint
	allocate( LU_x(npts,npts,2), LU_y(npts,npts,2), LU_z(npts,npts,2) )
	allocate( ipiv_x(npts,2), ipiv_y(npts,2), ipiv_z(npts,2) )

	!    x direction:
	if( nx > 1 ) then
		P = 2.d0 * Lx          ! periodicity length for even extensions
		do K=1,2
			if(K==1) then
				! near x=0
				a = 0.d0            ! shift for evaluating U_n(x-a)    here a=0
				istart = 0          ! array index-1  of first of npts x values to apply matching
			elseif(K==2) then
				! near x=Lx
				a = Lx              ! shift for evaluating U_n(x-a)    here a=Lx
				istart = nx - npts  ! array index-1  of first of npts x values to apply matching
			endif
		
			! fill the elements of the coefficient matrix 
			if(myid==0 .and. K==1) open(1,file='output/debug_data/Bcoeff_matrix_x_A')
			if(myid==0 .and. K==2) open(1,file='output/debug_data/Bcoeff_matrix_x_B')
			do i=1,npts
				xval = x(istart+i)
				do j=1,npts
					n = 2*j-1
					LU_x(i,j,K) = U_n(xval,a,P,n)
					if(myid==0) write(1,*) i,j,K,LU_x(i,j,K)
				enddo
			enddo
			if(myid==0) close(1)
		
			! decompose the coeff matrix into LU,  matrix is overwritten, pivot info in ipiv
			call DGETRF(npts,npts,LU_x(1,1,K),npts,ipiv_x(1,K),info)
			if( myid==0 .and. info .ne. 0 ) stop ' setup_factor_B_matrix problem, x direction '
		enddo
	endif
	
	!    y direction:
	if( ny > 1 ) then
		P = 2.d0 * Ly           ! periodicity length for even extensions
		do K=1,2
			if(K==1) then
				! near y=0
				a = 0.d0            ! shift for evaluating U_n(y-a)    here a=0
				istart = 0          ! array index-1  of first of npts y values to apply matching
			elseif(K==2) then
				! near y=Ly
				a = Ly              ! shift for evaluating U_n(y-a)    here a=Ly
				istart = ny - npts  ! array index-1  of first of npts y values to apply matching
			endif
		
			! fill the elements of the coefficient matrix
			if(myid==0 .and. K==1) open(1,file='output/debug_data/Bcoeff_matrix_y_A')
			if(myid==0 .and. K==2) open(1,file='output/debug_data/Bcoeff_matrix_y_B') 
			do i=1,npts
				yval = y(istart+i)
				do j=1,npts
					n = 2*j-1
					LU_y(i,j,K) = U_n(yval,a,P,n)
					if(myid==0) write(1,*) i,j,K,LU_y(i,j,K)
				enddo
			enddo
			if(myid==0) close(1)
		
			! decompose the coeff matrix into LU,  matrix is overwritten, pivot info in ipiv
			call DGETRF(npts,npts,LU_y(1,1,K),npts,ipiv_y(1,K),info)
			if( myid==0 .and. info .ne. 0 ) stop ' setup_factor_B_matrix problem, y direction '
		enddo	
	endif	
		
	!    z direction:
	if( nz > 1 ) then
		P = 2.d0 * Lz           ! periodicity length for even extensions
		do K=1,2
			if(K==1) then
				! near z=0
				a = 0.d0            ! shift for evaluating U_n(z-a)    here a=0
				istart = 0          ! array index-1  of first of npts z values to apply matching
			elseif(K==2) then
				! near z=Lz
				a = Lz              ! shift for evaluating U_n(z-a)    here a=Lz
				istart = nz - npts  ! array index-1  of first of npts z values to apply matching
			endif
		
			! fill the elements of the coefficient matrix
			if(myid==0 .and. K==1) open(1,file='output/debug_data/Bcoeff_matrix_z_A')
			if(myid==0 .and. K==2) open(1,file='output/debug_data/Bcoeff_matrix_z_B') 
			do i=1,npts
				zval = z(istart+i)
				do j=1,npts
					n = 2*j-1
					LU_z(i,j,K) = U_n(zval,a,P,n)
					if(myid==0) write(1,*) i,j,K,LU_z(i,j,K)
				enddo
			enddo
			if(myid==0) close(1)
		
			! decompose the coeff matrix into LU,  matrix is overwritten, pivot info in ipiv
			call DGETRF(npts,npts,LU_z(1,1,K),npts,ipiv_z(1,K),info)
			if( myid==0 .and. info .ne. 0 ) stop ' setup_factor_B_matrix problem, z direction '
		enddo
	endif
 return
END SUBROUTINE setup_factor_B_matrices
	
	
!-------------------------------------------------------------------------
SUBROUTINE solve_for_B_expansion_coeffs(LU,ipiv,f,nx,npts,x0,rhs,coeffs)
!--------------------------------------------------------------------------------------
!    fill the rhs vector and solve for the expansion coefficients for the U_n series
!    if the expansion point x0=0, then the coeffs are A1,A3,..A_Q and the series
!    should match f at the (Q+1)/2 points near x=0
!    if the expansion point x0=L, then the coeffs are B1,B3,..B_Q and the series
!    should match f at the (Q+1)/2 points near x=L
!    N.B. I'm passing in npts here instead of Q,   npts=(Q+1)/2
!--------------------------------------------------------------------------------------
	use mpi_params,           only: myid
	implicit none
	integer, intent(in)          :: nx,npts
	real(kind=8), intent(in)     :: lu(npts,npts),f(nx),x0
	integer, intent(in)          :: ipiv(npts)
	real(kind=8), intent(out)    :: coeffs(npts)
	real(kind=8), intent(inout)  :: rhs(npts)         ! input array is just work space
	integer                      :: istart, i
	integer                      :: nrhs=1, lda, ldb, info
	character(len=1)             :: trans='N'    ! solve Ax=b, not the transposed problem
	

	lda = npts
	ldb = npts
	
	if( x0==0.d0 ) then             ! expansion point at beginning of coord array
		istart = 0                  ! array index-1 of first x val to include
	else                            ! assume expansion point at end of coord array
		istart = nx - npts          ! array index-1 of first x val to include
	endif

	! fill the rhs vector
	do i=1,npts
		rhs(i) = f(istart+i)		
	enddo
	
	! solve for the npts expansion coefficients (rhs is overwritten with the soln)
	call DGETRS(trans, npts, nrhs, LU, lda, ipiv, rhs, ldb, info)	
	if( myid==0 .and. info .ne. 0 ) stop ' problem in solve_for_B_expansion_coeffs '
	coeffs(1:npts) = rhs(1:npts)

END SUBROUTINE solve_for_B_expansion_coeffs



!-------------------------------------------------------------------------
SUBROUTINE U_series_expansion(x,nx,npts,x0,coeffs,S,S_x)
!--------------------------------------------------------------------------------------
!    x is an nx vector of equally spaced points in [0,L]
!    construct the U_n series expansion with basis functions singular at x=x0
!    given x, Q the highest order in the (Q+1)/2 odd term expansion and 
!    the expansion coefficients coeffs[:]. Also construct the first derivative.
!    N.B. I'm passing in npts here instead of Q,   npts=(Q+1)/2
!--------------------------------------------------------------------------------------
	use mpi_params,           only: myid
	implicit none
	integer, intent(in)          :: nx,npts
	real(kind=8), intent(in)     :: x(nx),x0,coeffs(npts)
	real(kind=8), intent(inout)  :: S(nx),S_x(nx)
	real(kind=8)                 :: L,P
	integer                      :: i,j,n
	real(kind=8), external       :: U_n, ddx_U_n
	
	L = x(nx) - x(1)     ! domain length
	P = 2.d0*L           ! periodicity scale for even extensions
	S(:) = 0.d0
	S_x(:) = 0.d0
	
	do i=1,nx
		do j=1,npts
		 n = 2*j + 1
		 S(i)   = S(i)   + coeffs(j) *     U_n( x(i), x0, P, n )
		 S_x(i) = S_x(i) + coeffs(j) * ddx_U_n( x(i), x0, P, n )
		enddo
	enddo

 return
END SUBROUTINE U_series_expansion
	



!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION U_n(x,a,P,n)
!-------------------------------------------------------------------------
!    Evaluate the function U_n(x-a) defined in Eckhoff 98
!    MATHEMATICS OF COMPUTATION
!    Volume 67, Number 223, July 1998, pp. 1063-1087
!    The periodic in [0,1) function U_n(x) is a piecewise polynomial of degree n + 1
!     
!    a is shift, P is periodicity length for U_n(x), n is the order
!    N.B.  argument of the Bernoulli polynomials must be between 0 and 1
!
!      could c be optimized to reduce roundoff in the expansion sums????
!--------------------------------------------------------------------------------------
	implicit none
	integer, intent(in)       :: n
	real(kind=8), intent (in) :: x,a,P
	real(kind=8), external    :: B_poly
	integer, external         :: factorial     ! simple, for integers not too large
	real(kind=8)              :: xsi,c
	xsi = (x-a)/P
	c=-1.d0  !c = -(P)**n / factorial(n+1)        !!  I THINK NORMALIZATION OF U_n and ddx_U_n is arbitrary for my purposes
	!  wrap xsi as necessary for periodicity
	if( xsi < 0.d0 ) then 
		xsi = 1.d0 + xsi	
	endif
	
	! endpoint should have xsi=1 not xsi = 0 
	if( x==P/2.d0 .and. a==P/2.d0 ) xsi = 1.d0
	
	U_n = c*B_poly(xsi,n+1)					
 return
END FUNCTION U_n


!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION ddx_U_n(x,a,P,n)
!-------------------------------------------------------------------------
!    use the differentiation formula for the Bernoulli polynomials
!
!    Evaluate the function U_n(x-a) defined in Eckhoff 98
!    MATHEMATICS OF COMPUTATION
!    Volume 67, Number 223, July 1998, pp. 1063-1087
!    The periodic in [0,1) function U_n(x) is a piecewise polynomial of degree n + 1
!     
!    a is shift, L is periodicity length for U_n(x), n is the order
!    N.B.  argument of the Bernoulli polynomials must be between 0 and 1
!
!  U_n = c*B_n(xsi,n+1)                        U_n in terms of B_n+1	, xsi = (x-a)/L
!  d/dx U_n = c * d/dxsi B_n+1(xsi) * dxsi/dx
!  d/dxsi B_n+1(xsi) = (n+1)*B_n(xsi)
!--------------------------------------------------------------------------------------
	implicit none
	integer, intent(in)       :: n
	real(kind=8), intent (in) :: x,a,P
	real(kind=8), external    :: B_poly
	integer, external         :: factorial     ! simple, for integers not too large
	real(kind=8)              :: xsi,dxsi_dx,ddxsi,c
	xsi = (x-a)/P
	dxsi_dx = 1.d0/P       
	c=-1.d0  !c = -(P)**n / factorial(n+1)	      !!  I THINK NORMALIZATION OF U_n and ddx_U_n is arbitrary for my purposes
	!  wrap xsi as necessary for periodicity
	if( xsi < 0.d0 ) then 
		xsi = 1.d0 + xsi	
	endif
	
	! endpoint should have xsi=1 not xsi = 0 
	if( x==P/2.d0 .and. a==P/2.d0 ) xsi = 1.d0
	
							
	ddxsi = c * (n+1.d0)*B_poly(xsi,n)    ! d/dxsi of U_n(xsi)
	ddx_U_n = ddxsi * dxsi_dx             ! d/dx of U_n(x)
 return 
END FUNCTION ddx_U_n





!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION B_poly(x,n)
!-------------------------------------------------------------------------
!    Evaluate the nth Bernoulli polynomial B_n(x) using the 
!    "explicit formula" in Wikepedia https://en.wikipedia.org/wiki/Bernoulli_polynomials	  
!    ==> much more accurate than Komatsu and Pita Ruiz V. Math Commun. 21, 2016 pp. 127-140
!--------------------------------------------------------------------------------------
	implicit none
	integer, intent(in)      :: n
	real(kind=8),intent (in) :: x
	real(kind=8), external   :: bernoulli_number
	real(kind=8), external   :: bico
	integer                  :: k	
	
	B_poly = 0.d0
	do k=0,n	
		B_poly = B_poly + bico(n,k)*bernoulli_number(n-k)*x**k	
	enddo	
 return              ! returns B_n(x)
END FUNCTION B_poly


!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION ddx_B_n(x,n)
!-------------------------------------------------------------------------
!
!    Evaluate d/dx of the nth Bernoulli polynomial B_n(x) 
!-------------------------------------------------------------------------
	implicit none
	integer, intent(in)       :: n
	real(kind=8), intent (in) :: x
	real(kind=8), external    :: B_poly
	if( n > 0 ) then
		ddx_B_n = dfloat(n)*B_poly(x,n-1)    ! d/dx B_n(x) = n*B_n-1(x)
	else
		ddx_B_n = 0.d0
	endif				
 return 
END FUNCTION ddx_B_n


!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION mth_deriv_B_n(x,n,m)
!-------------------------------------------------------------------------------
!    Evaluate mth derivative of the nth Bernoulli polynomial B_n(x)
!    input params: n is Bernoulli polynomial index/order 
!                  m is the order of the derivative to evaluate at position x  
!-------------------------------------------------------------------------------
	implicit none
	integer, intent(in)       :: n,m
	real(kind=8), intent (in) :: x
	real(kind=8), external    :: B_poly
	integer, external         :: factorial     ! simple, for integers not too large	
	if( n >= m ) then
		mth_deriv_B_n = factorial(n)/factorial(n-m) * B_poly(x,n-m)   !  last term is B_(n-m)(x)
	else
		mth_deriv_B_n = 0.d0
	endif	
 return 
END FUNCTION mth_deriv_B_n








!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION bico(n,k)
!-------------------------------------------------------------------------
!
!   Compute the value "n choose k" = n!/( k!(n-k)! )  returns ans as float
!   uses factln   (from numerical recipes, special functions)
!-------------------------------------------------------------------------
	implicit none
	integer                :: k,n
	real(kind=8), external :: factln
	bico=nint(exp(factln(n)-factln(k)-factln(n-k)))
	!bico=exp(factln(n)-factln(k)-factln(n-k))
 return  !The nearest-integer function cleans up roundoff error for smaller values of n and k.
END FUNCTION bico 

!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION factln(n)
!-------------------------------------------------------------------------
!
!   Computes ln(n!)
!   uses gammln
!-------------------------------------------------------------------------
	implicit none
	integer                :: n
	real(kind=8), save     :: a(100)
	real(kind=8), external :: gammln
	data a/100*-1.d0/                            ! Initialize the table to negative values.
	
	if (n.lt.0) stop 'negative factorial in factln'
	if (n.le.99) then                              ! In range of the table.
		if (a(n+1).lt.0.) a(n+1)=gammln(n+1.d0)    ! If not already in the table, put it in.
		factln=a(n+1)
	else
		factln=gammln(n+1.d0)                      ! Out of range of the table.
	endif
	return
END FUNCTION factln


!-------------------------------------------------------------------------
REAL(kind=8) FUNCTION gammln(x)
!-------------------------------------------------------------------------
!
!   Compute log Gamma(x) for x > 0
!-------------------------------------------------------------------------
	integer            :: j
	real(kind=8)       :: ser,stp,tmp,x,y,cof(6)
	save               :: cof,stp
	data cof,stp/76.18009172947146d0,-86.50532032941677d0,                     \
                24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, \
                -.5395239384953d-5,2.5066282746310005d0/
	y=x
	tmp=x+5.5d0
	tmp=(x+0.5d0)*log(tmp)-tmp
	ser=1.000000000190015d0
	do j=1,6
		y=y+1.d0
		ser=ser+cof(j)/y
	enddo
	gammln=tmp+log(stp*ser/x)
	return
END FUNCTION gammln  
  
!-------------------------------------------------------------------------
REAL(Kind=8) FUNCTION bernoulli_number(n)
!-------------------------------------------------------------------------
!    Evaluate the nth Bernoulli number using the summation formula (2) in
!    Komatsu and Pita Ruiz V. Math Commun. 21, 2016 pp. 127-140
!-------------------------------------------------------------------------
	implicit none
	integer, intent(in)       :: n
	real(kind=8),external     :: bico
	integer                   :: i,j
	
	bernoulli_number = 0.d0
	do i=0,n
		do j=0,i
			bernoulli_number = bernoulli_number + (-1)**(i+j) * bico(n+1,j)/bico(n,i) * (i-j)**n
		enddo
	enddo
	bernoulli_number = bernoulli_number/(n+1.d0)
  	
 return
END FUNCTION bernoulli_number

!-------------------------------------------------------------------------
integer recursive function factorial(m) result (fac)
!-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m
    if( m < 0 ) stop 'argument of factorial function must be a positive integer'
    if( m == 0 ) then
      fac = 1
    else
      fac = m * factorial(m - 1)
    end if
  end function factorial


subroutine verify_deriv_BC
	use mpi_params,                       only: myid
	use differentiation_params,           only: Q
	use independent_variables,            only: x,y,z,Lx,Ly,Lz,nx,ny,nz,x_periodic,y_periodic,z_periodic
	implicit none 
	real(kind=8), allocatable                :: in(:),out(:)
	real(kind=8)                             :: pi, tol=1.e-8, diff
	integer                                  :: j, dir, nmin=64
	character(len=80), parameter             :: method='BC'
	logical                                  :: debug
	
	pi = 4.d0*atan(1.d0)	  

	!----------------------------------------------------
	! verify BC differentiation
	!----------------------------------------------------
	if( .not. x_periodic ) then
		allocate( in(nx),out(nx) )
		debug = .FALSE.
		in=0.d0 ; out=0.d0
 		do j=1,nx
 			in(j) = sin( pi*x(j)/Lx )   ! even extension has discontinuous derivs at x=0, x=Lx
 		enddo
 	
 		dir = 1
 		call deriv_BC(in,out,nx,dir,method,debug)  ! Bernoulli/Cosine derivative
 	
 		if(myid==0 .and. nx>nmin) then
 			do j=1,nx
 				diff = abs(out(j) - (pi/Lx)*cos( pi*x(j)/Lx ))
 				!write(0,*) x(j), in(j), out(j), (pi/Lx)*cos( pi*x(j)/Lx )
 				if( diff  > tol ) stop ' problem verifying x Bernoulli-cosine differentiation in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) '                         ....................... d/dx using Bernoulli-cosine method looks fine, tol=',tol
 		deallocate( in,out )
 	endif
 	
 	if( .not. y_periodic ) then
		allocate( in(ny),out(ny) )
		debug = .FALSE.
		in=0.d0 ; out=0.d0
 		do j=1,ny
 			in(j) = sin( pi*y(j)/Ly )   ! even extension has discontinuous derivs at y=0, y=Ly
 		enddo
 	
 		dir = 2
 		call deriv_BC(in,out,ny,dir,method,debug)  ! Bernoulli/Cosine derivative
 	
 		if(myid==0 .and. ny>nmin) then
 			do j=1,ny
 				diff = abs(out(j) - (pi/Ly)*cos( pi*y(j)/Ly ))
 				!write(0,*) y(j), in(j), out(j), (pi/Ly)*cos( pi*y(j)/Ly )
 				if( abs(diff) > tol ) stop ' problem verifying y Bernoulli-cosine differentiation in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) '                         ....................... d/dy using Bernoulli-cosine method looks fine, tol=',tol
 		deallocate( in,out )
 	endif
 	
 	if( .not. z_periodic ) then
 		allocate( in(nz),out(nz) )
 		debug = .TRUE.
		in=0.d0 ; out=0.d0
 		do j=1,nz
 			in(j) = sin( pi*z(j)/Lz )   ! even extension has discontinuous derivs at z=0, z=Lz
 		enddo
 	
 		dir = 3
 		call deriv_BC(in,out,nz,dir,method,debug)  ! Bernoulli/Cosine derivative
 	
 		if(myid==0 .and. nz>nmin) then
 			do j=1,nz
 				diff = abs(out(j) - (pi/Lz)*cos( pi*z(j)/Lz ))
 				!write(0,*) z(j), in(j), out(j), (pi/Lz)*cos( pi*z(j)/Lz ), diff
 				if( diff > tol ) stop ' problem verifying z Bernoulli-cosine differentiation in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) '                         ....................... d/dz using Bernoulli-cosine method looks fine, tol=',tol
 		deallocate( in,out )
 	endif 

 return
end subroutine verify_deriv_BC



