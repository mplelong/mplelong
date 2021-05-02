!-----------------------------------------------------------------------------------------
!
!             solver utilities.... routines that know about and operate on
!             data arranged in XBLOCK, YBLOCK and ZBLOCK storage schemes
!             or some other information specific to flow_solve
!
!-----------------------------------------------------------------------------------------

subroutine ddx(f,df,order,Qval)
	!-------------------------------------------------------------------------
	! assume f and df are decomposed in XBLOCK decomposition
	!-------------------------------------------------------------------------
	use mpi_params,             only: myid
 	use independent_variables,  only: nx
 	use decomposition_params,   only: array_size,IDIM,JDIM,KDIM,XBLOCK 
 	implicit none 
 	integer                        ::  j,k,dir=1
 	real(kind=8), intent(in)       ::  f( array_size(IDIM,XBLOCK,myid),    &
                                          array_size(JDIM,XBLOCK,myid),    &
                                          array_size(KDIM,XBLOCK,myid)  )
	real(kind=8), intent(inout)    :: df( array_size(IDIM,XBLOCK,myid),    &
                                          array_size(JDIM,XBLOCK,myid),    &
                                          array_size(KDIM,XBLOCK,myid)  )
 	integer, intent(in)            :: order,Qval
 	logical,save                   :: first_entry=.TRUE. , debug=.FALSE.
	!--------------------------------------------------------------------------

	if(first_entry) then
 		if( nx .ne. array_size(IDIM,XBLOCK,myid) ) stop 'XBLOCK decomposition error, ddx'
 		first_entry=.FALSE.
	endif

	if(nx==1) then
 		df(:,:,:) = 0.d0
 		return
	endif

	!--------------------------------------------------------------------------
	! loop through 2nd & 3rd array indices, perform stride 1
	! global differentiation operation, with data local to myid
	!--------------------------------------------------------------------------  
 	do k=1,array_size(KDIM,XBLOCK,myid)  
  		do j=1,array_size(JDIM,XBLOCK,myid)   		
  			call deriv_BC(f(1,j,k),df(1,j,k),nx,dir,Qval,debug)  ! Bernoulli/Cosine derivative
  		enddo
 	enddo
 
 return
end subroutine ddx


subroutine ddy(f,df,order,Qval)
	!-------------------------------------------------------------------------
	! assume f and df are decomposed in YBLOCK decomposition
	!-------------------------------------------------------------------------
	use mpi_params,             only: myid
	use independent_variables,  only: ny
 	use decomposition_params,   only: array_size,IDIM,JDIM,KDIM,YBLOCK
 	implicit none 
 	integer                        :: i,k,dir=2
 	real(kind=8), intent(in)       ::  f( array_size(IDIM,YBLOCK,myid),    &
                                          array_size(JDIM,YBLOCK,myid),    &
                                          array_size(KDIM,YBLOCK,myid)  )
 	real(kind=8), intent(out)      :: df( array_size(IDIM,YBLOCK,myid),    &
                                          array_size(JDIM,YBLOCK,myid),    &
                                          array_size(KDIM,YBLOCK,myid)  )
 	integer, intent(in)            :: order,Qval 
 	logical,save                   :: first_entry=.TRUE. , debug=.FALSE.
	!--------------------------------------------------------------------------
	if(first_entry) then
		if( ny .ne. array_size(IDIM,YBLOCK,myid) ) stop 'YBLOCK decomposition error, ddy'
 		first_entry=.FALSE.
	endif

	if(ny==1) then
 		df(:,:,:)=0.d0
 		return
	endif
 	
	!--------------------------------------------------------------------------
	! loop through x and z array indices, perform stride 1
	! global differentiation operation with data local to myid
	!--------------------------------------------------------------------------
 	do k=1,array_size(KDIM,YBLOCK,myid)  
  		do i=1,array_size(JDIM,YBLOCK,myid)
  			call deriv_BC(f(1,i,k),df(1,i,k),ny,dir,Qval,debug)  ! Bernoulli/Cosine derivative 
 		enddo
 	enddo
 
 return
end subroutine ddy


subroutine ddz(f,df,order,Qval)
	!-------------------------------------------------------------------------
	! assume f and df are decomposed in ZBLOCK decomposition (z,x,y)
	!-------------------------------------------------------------------------
	use mpi_params,             only: myid,comm,ierr
 	use independent_variables,  only: nz
 	use decomposition_params,   only: array_size,IDIM,JDIM,KDIM,ZBLOCK 
 	implicit none 
 	integer                        :: j,k,dir=3
 	real(kind=8)                   ::  f( array_size(IDIM,ZBLOCK,myid),    &
                                          array_size(JDIM,ZBLOCK,myid),    &
                                          array_size(KDIM,ZBLOCK,myid)  )
 	real(kind=8)                   :: df( array_size(IDIM,ZBLOCK,myid),    &
                                          array_size(JDIM,ZBLOCK,myid),    &
                                          array_size(KDIM,ZBLOCK,myid)  )
 	integer, intent(in)            :: order,Qval
 	logical,save                   :: first_entry=.TRUE. , debug=.FALSE.
	!--------------------------------------------------------------------------

	if(first_entry) then
 		if( nz .ne. array_size(IDIM,ZBLOCK,myid) ) stop 'ZBLOCK decomposition error, ddz'
 		first_entry=.FALSE.
	endif

	if(nz==1) then
 		df(:,:,:)=0.d0
 		return
	endif

	!--------------------------------------------------------------------------
	! loop through 2nd & 3rd array indices, perform stride 1
	! global differentiation operation, with data local to myid
	!--------------------------------------------------------------------------
 	do k=1,array_size(KDIM,ZBLOCK,myid)  
  		do j=1,array_size(JDIM,ZBLOCK,myid)  
  			call deriv_BC(f(1,j,k),df(1,j,k),nz,dir,Qval,debug)  ! Bernoulli/Cosine derivative
  		enddo
 	enddo

 return
end subroutine ddz


subroutine divergence(u,v,w,div)
	!------------------------------------------------------------
	!  Compute divergence of [u,v,w] using BernoulliCosine method
	!  Input and output data arrays arranged in YBLOCK format.
	!  Uses tmpY(:,:,:,6) for work space.
	!------------------------------------------------------------
	use mpi_params,             only: myid,comm,ierr
 	use decomposition_params
 	use independent_variables,  only: nx,ny,nz,x_periodic,y_periodic,z_periodic
 	use intermediate_variables, only: tmpX,tmpY,tmpZ
 	use differentiation_params, only: Q            ! (Q+1)/2 terms in Bernoulli series
 
 	implicit none
 	include 'mpif.h'
 	integer,parameter             ::   order=1
 	integer                       ::   Qval
 	real(kind=8), intent(in)      ::   u( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )
                                       
 	real(kind=8), intent(in)      ::   v( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )
                                       
 	real(kind=8), intent(in)      ::   w( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )
                                       
 	real(kind=8), intent(out)     :: div( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  ) 

 	!----------------------------------------
 	!    div <--- dv/dy
 	!----------------------------------------
 	if( y_periodic ) then
 		Qval = 0
 	else	
 		Qval = Q
 	endif
 	
	call ddy( v, div, order, Qval )                        ! ==> dv/dy in YBLOCK format
	
 
 	!----------------------------------------
 	!    div <--- dv/dy + dw/dz
 	!----------------------------------------
 	if( z_periodic ) then
 		Qval = 0
 	else	
 		Qval = Q
 	endif
 	
  	call yblock_2_zblock(w(1,1,1),tmpZ(1,1,1,1))   	
  	call ddz( tmpZ(1,1,1,1), tmpZ(1,1,1,2), order, Qval )   ! ==> dw/dz in ZBLOCK format  	
  	call zblock_2_yblock(tmpZ(1,1,1,2),tmpY(1,1,1,6))       ! ==> dw/dz in YBLOCK format
  	
  	div(:,:,:) = div(:,:,:) + tmpY(:,:,:,6)
 	
 
 	!----------------------------------------
 	!    div <--- dv/dy + dw/dz + dudx
 	!----------------------------------------
 	if( x_periodic ) then
 		Qval = 0
 	else	
 		Qval = Q
 	endif
 	
  	call yblock_2_xblock(u,tmpX) 		
  	call ddx( tmpX, tmpX(1,1,1,2), order, Qval )           ! ==> du/dx in XBLOCK format 		
  	call xblock_2_yblock(tmpX(1,1,1,2),tmpY(1,1,1,6))      ! ==> du/dx in YBLOCK format
  	
  	div(:,:,:) = div(:,:,:) + tmpY(:,:,:,6)
  
 return
end subroutine divergence

subroutine test_divergence
!---------------------------------------------------------------
! Simple routine to test the divergence routine
!  to pass the test, the transpose and the differentiation
!  routines must work. This routine should be called from
!  preliminary tasks.
!---------------------------------------------------------------
	use mpi_params,              only: myid,comm,ierr
	use independent_variables,   only: Lx,Ly,Lz,nx,nz
	use intermediate_variables,  only: tmpY,tmpZ
	use decomposition_params
	use etc
	implicit none
	real(kind=8), allocatable      :: x(:),y(:),z(:) ! local portions
	integer                        :: i,j,k,locnx,ny,locnz,nmin=64
	real(kind=8)                   :: ans,pi,kval,diff,tol=1.d-9
	
	
	locnx = array_size(JDIM,YBLOCK,myid)
	ny    = array_size(IDIM,YBLOCK,myid)
	locnz = array_size(KDIM,YBLOCK,myid)
	
	if( nx>1 .AND. ny>1 .AND. nz>1 ) then
		! don't do a 3d test with strict tolerance for small problems
		if( nx< nmin .OR. ny < nmin .OR. nz < nmin ) return
	endif
		
	!--------------------
	! local x,y,z vals
	!--------------------
	allocate( x(locnx), y(ny), z(locnz) )
	call get_my_xvals(x,YBLOCK,myid)
	call get_my_yvals(y,YBLOCK,myid)
	call get_my_zvals(z,YBLOCK,myid)		
	
	!---------------------------------------------------------------------------
	! define a test vector:  
	!                          [ (x/Lx)^2, cos(kval*y) + (y/Ly)^2, (z/Lz)^2 ]
	!    ==>  div = 2( x/Lx^2 + y/Ly^2 + z/Lz^2 ) - kval*sin(kval*y)
	!---------------------------------------------------------------------------
	pi = 4.d0*atan(1.d0)
	kval = 0.d0 * (2.*pi/Ly)
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				tmpY(j,i,k,1) = (x(i)/Lx)**2
				tmpY(j,i,k,2) = (y(j)/Ly)**2 + cos( kval*y(j) )
				tmpY(j,i,k,3) = (z(k)/Lz)**2
			enddo
		enddo
	enddo
	
	!---------------------------------------------------------------------------
	!  compute the divergence, store result in array 4
	!---------------------------------------------------------------------------
	call divergence(tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),tmpY(1,1,1,4))
	
	!---------------------------------------------------------------------------
	!  compare computed divergence to exact values
	!---------------------------------------------------------------------------
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				ans = 2.d0*( x(i)/Lx**2 + y(j)/Ly**2 + z(k)/Lz**2 ) - kval*sin(kval*y(j))
				diff = ans - tmpY(j,i,k,4)
				!if(i==1 .and. j==1) write(0,*) z(k)/Lz, ans, tmpY(j,i,k,4), diff				
				if( abs( diff ) > tol ) stop '    test of divergence routine failed '
			enddo
		enddo
	enddo
	call mpi_barrier(comm,ierr) 
	
	if(myid==0) then
  		message = '...  Analytical test of divergence routine successful ... tol=1.d-9'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 return
end subroutine test_divergence


subroutine gradient(f,fx,fy,fz,Q0)
	!------------------------------------------------------------
	!  Compute spatial gradient of f using Bernoulli/Cos method.
	!  Input and output data arrays arranged in YBLOCK format.
	!  Arg for Qval allows use of standard Cos gradient for phi.
	!------------------------------------------------------------
	use mpi_params,             only: myid,comm,ierr
	use decomposition_params
	use intermediate_variables
	use independent_variables
	implicit none
	integer                       ::   order=1
	integer                       ::   Q0,Qval
	real(kind=8)                  ::   f( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )
                                       
	real(kind=8)                  ::  fx( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )
                                       
	real(kind=8)                  ::  fy( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )
                                       
	real(kind=8)                  ::  fz( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid)  )    

 	if( array_size(IDIM,YBLOCK,myid) > 1 ) then
 		if( y_periodic ) then
 			Qval = 0
 		else	
 			Qval = Q0
 		endif
 	
  		call ddy( f, fy, order, Qval )
 	else
  		fy = 0.d0
 	endif
 	

 	if( nz > 1 ) then
 		if( z_periodic ) then
 			Qval = 0
 		else	
 			Qval = Q0
 		endif
 		
  		call yblock_2_zblock(f,tmpZ)  
  		call ddz( tmpZ, tmpZ(1,1,1,2), order, Qval )
  		call zblock_2_yblock(tmpZ(1,1,1,2),fz)
 	else
  		fz = 0.d0
 	endif


 	if( nx > 1 ) then
 		if( x_periodic ) then
 			Qval = 0
 		else	
 			Qval = Q0
 		endif
 		
  		call yblock_2_xblock(f,tmpX)
  		call ddx( tmpX, tmpX(1,1,1,2), order, Qval )
  		call xblock_2_yblock(tmpX(1,1,1,2),fx)
 	else
  		fx = 0.d0
 	endif
 
 return
end subroutine gradient

subroutine test_gradient
	!---------------------------------------------------------------
	! Simple routine to test the gradient routine
	!  to pass the test, the transpose and the differentiation
	!  routines must work. This routine should be called from
	!  preliminary tasks.
	!---------------------------------------------------------------
	use mpi_params,              only: myid,comm,ierr
	use independent_variables,   only: Lx,Ly,Lz,nx,nz
	use intermediate_variables,  only: tmpY,tmpX,tmpZ
	use decomposition_params
	use differentiation_params,  only: Q
	use etc
	implicit none
	integer                        :: i,j,k,nmin=64
	real(kind=8),allocatable,save  :: x(:),y(:),z(:) ! local portions
	integer, save                  :: locnx,ny,locnz
	real(kind=8)                   :: ans,diff,tol=1.d-8
	logical, save                  :: first_entry=.TRUE.
	
	if( first_entry) then
		locnx = array_size(JDIM,YBLOCK,myid)
		ny    = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		
		!--------------------
		! local x,y,z vals
		!--------------------
		allocate( x(locnx), y(ny), z(locnz) )
		x = 0.d0 ; y=0.d0; z=0.d0
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)
		 			
		first_entry = .FALSE.
	endif
	
	if( nx>1 .AND. ny>1 .AND. nz>1 ) then
		! don't do a 3d test with strict tolerance for small problems
		if( nx< nmin .OR. ny < nmin .OR. nz < nmin ) return
	endif
	
	!--------------------------------------------------------------------
	! define a test function:  f = (x/Lx)^2 + (y/Ly)^2 + (z/Lz)^2
	!                          f_x = 2x/Lx^2  f_y=2y/Ly^2   f_z=2z/Lz^2
	!--------------------------------------------------------------------
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				tmpY(j,i,k,4) = (x(i)/Lx)**2 + (y(j)/Ly)**2 + (z(k)/Lz)**2
			enddo
		enddo
	enddo
	
	!-------------------------------------------------------------------------
	!  compute the gradient, put components of the answer in tmpY(:,:,:,1-3)
	!-------------------------------------------------------------------------
	call gradient(tmpY(1,1,1,4),tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),Q)
	
	!-------------------------------------------------------------------------
	!  compare computed gradient to exact values
	!-------------------------------------------------------------------------
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				ans = 2.d0*( x(i)/Lx**2  )
				diff = abs(ans - tmpY(j,i,k,1))
				if( diff > tol ) stop '    test of gradient routine failed, x direction '
				
				ans = 2.d0*( y(j)/Ly**2  )
				diff = abs(ans - tmpY(j,i,k,2))
				if( diff > tol ) stop '    test of gradient routine failed, y direction '
				
				ans = 2.d0*( z(k)/Lz**2  )
				diff = abs(ans - tmpY(j,i,k,3))
				!if(i==1 .and. j==1) write(0,*) ans,tmpY(j,i,k,3),diff
				if( diff > tol ) stop '    test of gradient routine failed, z direction '
				
			enddo
		enddo
	enddo
	call mpi_barrier(comm,ierr)
	
	if(myid==0) then
  		message = '...  Analytical test of gradient routine successful   ...  tol=1.d-10'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
	
 return
end subroutine test_gradient



subroutine udotgradf(u,v,w,f,ans)
	!---------------------------------------------------
	! Assume input and output arrays in YBLOCK format
	!   uses tmpY(:,:,:,1-3) to hold f_x, f_y and f_z
	!---------------------------------------------------
	use mpi_params,             only: myid
 	use decomposition_params
 	use intermediate_variables, only: tmpY
 	use differentiation_params, only: Q
 	implicit none
 	logical                   :: debug=.TRUE.
	real(kind=8)              :: ans(array_size(IDIM,YBLOCK,myid),     &
                                     array_size(JDIM,YBLOCK,myid),     &
                                     array_size(KDIM,YBLOCK,myid) )

 	real(kind=8)              ::   u(array_size(IDIM,YBLOCK,myid),     &
                                     array_size(JDIM,YBLOCK,myid),     &
                                     array_size(KDIM,YBLOCK,myid) )

 	real(kind=8)              ::   v(array_size(IDIM,YBLOCK,myid),     &
                                     array_size(JDIM,YBLOCK,myid),     &
                                     array_size(KDIM,YBLOCK,myid) )

 	real(kind=8)              ::   w(array_size(IDIM,YBLOCK,myid),     &
                                     array_size(JDIM,YBLOCK,myid),     &
                                     array_size(KDIM,YBLOCK,myid) )

 	real(kind=8)              ::   f(array_size(IDIM,YBLOCK,myid),     &
                                     array_size(JDIM,YBLOCK,myid),     &
                                     array_size(KDIM,YBLOCK,myid) )

	!---------------------------------------------------
	! compute grad f, store results in tmpY 1,2 and 3
	!---------------------------------------------------
 	call gradient(f,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),Q)

	!---------------------------------------------------
	! do the dot product
	!---------------------------------------------------
 	ans(:,:,:) = u(:,:,:)*tmpY(:,:,:,1) + v(:,:,:)*tmpY(:,:,:,2) + w(:,:,:)*tmpY(:,:,:,3)

 return
end subroutine udotgradf

subroutine test_udotgradf
!---------------------------------------------------------------
! Simple routine to test the udotgradf routine
!  to pass the test, the transpose and the differentiation
!  routines must work. This routine should be called from
!  preliminary tasks.
!---------------------------------------------------------------
	use mpi_params,              only: myid,comm,ierr
	use independent_variables,   only: Lx,Ly,Lz,nx,nz
	use intermediate_variables,  only: tmpY
	use decomposition_params
	use differentiation_params,  only: Q
	use etc
	implicit none
	real(kind=8), allocatable      :: x(:),y(:),z(:) ! local portions
	integer                        :: i,j,k,locnx,ny,locnz,nmin=64
	real(kind=8)                   :: ans,diff,tol=1.d-8
	logical, save                  :: first_entry=.TRUE.
	
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
		first_entry = .FALSE.
	endif
	
	if( nx>1 .AND. ny>1 .AND. nz>1 ) then
		! don't do a 3d test with strict tolerance for small problems
		if( nx< nmin .OR. ny < nmin .OR. nz < nmin ) return
	endif
	
	!----------------------------------------------------------------------
	! define a test function:  f = (x/Lx)^2 + (y/Ly)^2 + (z/Lz)^2
	!                          f_x = 2x/Lx^2  f_y=2y/Ly^2   f_z=2z/Lz^2
	!----------------------------------------------------------------------
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				tmpY(j,i,k,4) = (x(i)/Lx)**2 + (y(j)/Ly)**2 + (z(k)/Lz)**2  ! f
				tmpY(j,i,k,5) = 3.d0                                        ! u = v = w
			enddo
		enddo
	enddo
	
	!----------------------------------------------------------------
	! udotgradf uses tmpY(:,:,:,1-3) to store f_x, f_y and f_z
	! after computing them, here to avoid these arrays use u=v=w=3.d0,  
	! compute u dot grad f and store answer in tmpY(1,1,1,6)
	!----------------------------------------------------------------
	call udotgradf(tmpY(1,1,1,5),tmpY(1,1,1,5),tmpY(1,1,1,5),tmpY(1,1,1,4),tmpY(1,1,1,6))
		
	!  compare computed to exact values
	do k=1,locnz
		do i=1,locnx
			do j=1,ny
				ans = 3.d0*( 2.d0*( x(i)/Lx**2  ) ) &
				    + 3.d0*( 2.d0*( y(j)/Ly**2  ) ) & 
				    + 3.d0*( 2.d0*( z(k)/Lz**2  ) )		
				diff = abs(ans - tmpY(j,i,k,6))	
				if(j==0 .and. k==0 ) write(0,*) j,i,k,ans,tmpY(j,i,k,6)	    
				if( abs( ans - tmpY(j,i,k,6) ) > tol ) stop '    test of udotgradf routine failed with tol=1.d-10'
			enddo
		enddo
	enddo
	
	if(myid==0) then
  		message = '...  Analytical test of udotgradf routine successful ...  tol=1.d-10'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
	
 return
end subroutine test_udotgradf


subroutine divgradf(f,g2f)
	use mpi_params,             only: myid
 	use decomposition_params
	use intermediate_variables, only: tmpY
 	use differentiation_params, only: Q
 	implicit none
 	real(kind=8)              ::   f(array_size(IDIM,YBLOCK,myid),     &
                                     array_size(JDIM,YBLOCK,myid),     &
                                     array_size(KDIM,YBLOCK,myid) )
    real(kind=8)              :: g2f(array_size(IDIM,YBLOCK,myid),     &
                                     array_size(JDIM,YBLOCK,myid),     &
                                     array_size(KDIM,YBLOCK,myid) )

	! compute gradient of f, store in tmpY 1,2 & 3
	call gradient(f,tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),Q) 
	
	! compute div of grad f 
	call divergence(tmpY(1,1,1,1),tmpY(1,1,1,2),tmpY(1,1,1,3),g2f) 	
end subroutine divgradf



subroutine initialize_fourier_stuff
	!------------------------------------------------------------------
	!  Allocate and fill wavenumber and filter arrays
	!  Create the required FFTW3 plans
	!  Store results in differentiation_params module
	!------------------------------------------------------------------
	use mpi_params,                 only: myid,comm,ierr
	use independent_variables
	use differentiation_params
	use fourier_differentiation_tools 
	implicit none
	integer                            :: i,j,k
	character(len=80)                  :: exp_type
	integer                            :: dim,dir,order
	integer(kind=8)                    :: plan_i,plans(2)
	real(kind=8)                       :: tol = 1.d-8, pi
	real(kind=8), allocatable          :: in(:), out(:)

	allocate( kx(nx),ky(ny),kz(nz) )
	kx=0.d0 ; ky=0.d0; kz=0.d0
 	allocate( kxfilter(nx), kyfilter(ny), kzfilter(nz) )
 	kxfilter=0.d0 ; kyfilter=0.d0 ; kzfilter=0.d0
 	
 	dim = 1   ! x
 	if( x_periodic ) then
 		exp_type = 'fourier'
		call fourier_init(exp_type,nx,Lx,kx,kxfilter,fourier_plan(dim,1),fourier_plan(dim,2))
		fourier_done(dim) = .TRUE.
 	else
 		exp_type = 'cos'
		call fourier_init(exp_type,nx,Lx,kx,kxfilter,cos_plan(dim),plan_i)
		exp_type = 'sin'
		call fourier_init(exp_type,nx,Lx,kx,kxfilter,sin_plan(dim),plan_i)
		cos_done(dim) = .TRUE.
		sin_done(dim) = .TRUE.
	endif
	if(myid==0) then
		open(1,file='output/debug_data/x_wavenumbers')
		do i=1,nx
			write(1,*) i,kx(i),kxfilter(i)
		enddo
		close(1)
	endif
	
	dim = 2   ! y
	if( y_periodic ) then
		exp_type = 'fourier'
		call fourier_init(exp_type,ny,Ly,ky,kyfilter,fourier_plan(dim,1),fourier_plan(dim,2))
		fourier_done(dim) = .TRUE.
	else
 		exp_type = 'cos'
		call fourier_init(exp_type,ny,Ly,ky,kyfilter,cos_plan(dim),plan_i)
		exp_type = 'sin'
		call fourier_init(exp_type,ny,Ly,ky,kyfilter,sin_plan(dim),plan_i)
		cos_done(dim) = .TRUE.
		sin_done(dim) = .TRUE.
	endif
	if(myid==0) then
		open(1,file='output/debug_data/y_wavenumbers')
		do j=1,ny
			write(1,*) j,ky(j),kyfilter(j)
		enddo
		close(1)
	endif
	
	dim = 3   ! z
	if( z_periodic ) then
		exp_type = 'fourier'
		call fourier_init(exp_type,nz,Lz,kz,kzfilter,fourier_plan(dim,1),fourier_plan(dim,2))
		fourier_done(dim) = .TRUE.
	else
 		exp_type = 'cos'
		call fourier_init(exp_type,nz,Lz,kz,kzfilter,cos_plan(dim),plan_i)
		exp_type = 'sin'
		call fourier_init(exp_type,nz,Lz,kz,kzfilter,sin_plan(dim),plan_i)
		cos_done(dim) = .TRUE.
		sin_done(dim) = .TRUE.
	endif
	if(myid==0) then
		open(1,file='output/debug_data/z_wavenumbers')
		do k=1,nz
			write(1,*) k,kz(k),kzfilter(k)
		enddo
		close(1)
	endif
	
	
	!----------------------------------------------------
	! verify fourier or cosine differentiation
	!----------------------------------------------------
	pi = 4.d0*atan(1.d0)
	
	allocate( in(nx),out(nx) )
		in=0.d0 ; out=0.d0
 		do j=1,nx
 			in(j) = cos( 2.*pi*x(j)/Lx )
 		enddo
 	
 		dir = 1
 		if( x_periodic ) then
 			plans(1) = fourier_plan(dir,1)
 			plans(2) = fourier_plan(dir,2)
 			exp_type = 'fourier'
 		else
 			plans(1) = cos_plan(dir)
 			plans(2) = sin_plan(dir)
 			exp_type = 'cos'
 		endif
 		order = 1
 		!call fourier_deriv(in,out,ny,order,exp_type,ky,kyfilter,tmpY,plans)  ! low level call
 		call differentiate_fcs(in,out,nx,dir,exp_type,order)                  ! higher level call
 	 	 	
 		if(myid==0) then
 			do j=1,nx
 				!write(0,*) x(j), in(j), out(j), (-2.*pi/Lx)*sin( 2.*pi*x(j)/Lx )
 				if( abs(out(j) - (-2.*pi/Lx)*sin( 2.*pi*x(j)/Lx )) > tol ) stop ' problem verifying cosine differentiation in x in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) ' ....................... d/dx using term by term differentiation of cosine series looks fine, tol=',tol
 	deallocate( in,out )
 	
	allocate( in(ny),out(ny) )
		in=0.d0 ; out=0.d0
 		do j=1,ny
 			in(j) = cos( 2.*pi*y(j)/Ly )
 		enddo
 	
 		dir = 2
 		if( y_periodic ) then
 			plans(1) = fourier_plan(dir,1)
 			plans(2) = fourier_plan(dir,2)
 			exp_type = 'fourier'
 		else
 			plans(1) = cos_plan(dir)
 			plans(2) = sin_plan(dir)
 			exp_type = 'cos'
 		endif
 		order = 1
 		!call fourier_deriv(in,out,ny,order,exp_type,ky,kyfilter,tmpY,plans)  ! low level call
 		call differentiate_fcs(in,out,ny,dir,exp_type,order)                  ! higher level call
 	 	 	
 		if(myid==0) then
 			do j=1,ny
 				!write(0,*) y(j), in(j), out(j), (-2.*pi/Ly)*sin( 2.*pi*y(j)/Ly )
 				if( abs(out(j) - (-2.*pi/Ly)*sin( 2.*pi*y(j)/Ly )) > tol ) stop ' problem verifying cosine differentiation in y in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) ' ....................... d/dy using term by term differentiation of cosine series looks fine, tol=',tol
 	deallocate( in,out )
 	
 	allocate( in(nz),out(nz) )
		in=0.d0 ; out=0.d0
 		do j=1,nz
 			in(j) = cos( 2.*pi*z(j)/Lz )
 		enddo
 	
 		dir = 3
 		if( z_periodic ) then
 			plans(1) = fourier_plan(dir,1)
 			plans(2) = fourier_plan(dir,2)
 			exp_type = 'fourier'
 		else
 			plans(1) = cos_plan(dir)
 			plans(2) = sin_plan(dir)
 			exp_type = 'cos'
 		endif
 		order = 1
 		!call fourier_deriv(in,out,ny,order,exp_type,ky,kyfilter,tmpY,plans)  ! low level call
 		call differentiate_fcs(in,out,nz,dir,exp_type,order)                  ! higher level call
 		 	 	 	
 		if(myid==0) then
 			do j=1,nz
 				!write(0,*) y(j), in(j), out(j), (-2.*pi/Ly)*sin( 2.*pi*y(j)/Ly )
 				if( abs(out(j) - (-2.*pi/Lz)*sin( 2.*pi*z(j)/Lz )) > tol ) stop ' problem verifying cosine differentiation in z in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) ' ....................... d/dz using term by term differentiation of cosine series looks fine, tol=',tol
 	deallocate( in,out )
		
 return
end subroutine initialize_fourier_stuff




subroutine cfl
	use etc
	use mpi_params
	use independent_variables,  only: dt,x,y,z,nx,ny,nz,t_secs
	use dependent_variables
	implicit none
	include 'mpif.h'
	real(kind=8),save :: dx,dy,dz
	real(kind=8)      :: cfl_x,cfl_y,cfl_z
	real(kind=8)      :: umax,vmax,wmax,global_max
	integer           :: icount,root_pid
	logical           :: shutdown=.FALSE.
	logical, save     :: first_entry=.TRUE.	
	if( first_entry ) then
		if( nx > 1 ) then
			dx = x(2)-x(1)
		else
			dx = 1.e16     ! essentially infinite	
		endif
		if( ny > 1 ) then
			dy = y(2)-y(1)
		else
			dy = 1.e16     ! essentially infinite	
		endif
		if( nz > 1 ) then
			dz = z(2)-z(1)
		else
			dz = 1.e16     ! essentially infinite	
		endif
		first_entry = .FALSE.	
	endif
	icount=1         ! compute max of single value
	root_pid=0       ! to return mpi_reduce results to pid=0

 	!------------------------------
 	! local value of max speeds
 	!------------------------------
  	umax = MAXVAL( ABS(u) )
  	vmax = MAXVAL( ABS(v) )
  	wmax = MAXVAL( ABS(w) )

 	!------------------------------
 	! local max CFL values
 	!------------------------------
 	cfl_x = umax*dt/dx
 	cfl_y = vmax*dt/dy
 	cfl_z = wmax*dt/dz
 	
 	!----------------------------------------------------
 	!  find the maximum values across all the processors
 	!  use ALLreduce so each processor gets the max value
 	!----------------------------------------------------
 	call MPI_ALLREDUCE(cfl_x,global_max,icount,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
 	cfl_x = global_max

 	call MPI_ALLREDUCE(cfl_y,global_max,icount,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
 	cfl_y = global_max

 	call MPI_ALLREDUCE(cfl_z,global_max,icount,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
 	cfl_z = global_max

 	!-------------------------------------------------------
 	! only processor 0 writes the result to the screen log
 	!-------------------------------------------------------
 	if( myid==0 ) then
  		write(0,*) '...............     ADVECTIVE CFL ratios '
		write(0,*) '...............              x :  ',cfl_x
  		write(0,*) '...............              y :  ',cfl_y
  		write(0,*) '...............              z :  ',cfl_z
  		open(1,file='output/cfl.dat',position='append')
   			write(1,1001) t_secs,cfl_x,cfl_y,cfl_z
  		close(1)
 	endif
1001 format (4e18.8)

	!------------------------------------------------------------------------------
	! any and all processors execute fortran "STOP" if there is a cfl problem
	! when only root processor issued stop, there were hangs on some machines...
	!------------------------------------------------------------------------------
  	if( cfl_x > 1 ) stop 'CFL_x > 1 incipient instability '
  	if( cfl_y > 1 ) stop 'CFL_y > 1 incipient instability '
  	if( cfl_z > 1 ) stop 'CFL_z > 1 incipient instability '
  	if( istep > 0 .and. cfl_x==0. .and. cfl_y==0. .and. cfl_z==0. )  then
   		write(0,*) 'all CFL numbers zero: probably post-instability myid:',myid
   		shutdown = .TRUE.
   	endif
   	call mpi_barrier(comm,ierr)
   	if( shutdown ) call quit_mpi
 return
end subroutine cfl



subroutine toggle
	use etc
	use independent_variables,  only: t_secs,dt,t0,tf,tn,tnp1
	use mpi_params,             only: myid
	use methods_params,         only: AB_order
	use timing
	implicit none
	integer                        :: itmp
	integer,parameter              :: iprint=25
	logical, save                  :: first_entry=.TRUE.

 	if( myid .eq. 0 .and. mod(istep,iprint)==0 ) then
  		write(0,*) '..........................................................'
  		write(0,*) '............................     Toggling, istep  ',istep
  		write(0,*) '..........................................................'
 	endif
 	
 	if( first_entry ) then
 		iend = (tf-t0)/dt + 1
 		first_entry=.FALSE.
 	endif

 	t_secs = t_secs + dt
 	tn = t_secs
 	tnp1 = tn + dt      ! used to evaluate bvals at next time step

	! Toggle M cycle pointers for AB rhs fields
 	if( AB_order .eq. 4 ) then
		itmp=MM3
		MM3=MM2
		MM2=MM1
		MM1=MM0
		MM0=itmp
	elseif( AB_order .eq. 3 ) then
		itmp=MM2
		MM2=MM1
		MM1=MM0
		MM0=itmp
	elseif( AB_order .eq. 2 ) then
		itmp=MM1
		MM1=MM0
		MM0=itmp
	endif

	t_end_time_step = mpi_wtime()
	t_time_step = t_end_time_step - t_start_time_step
	t_total = t_total + t_time_step

	if( mod(istep,iprint)==0 ) call cfl

	istep=istep+1
	if(istep > iend) then
		if(myid==0) then
			write(0,*) ' time for last time step:   ',t_time_step
			write(0,*) '              total time:   ',t_total
		endif
		call finalize
	endif
 return
end subroutine toggle

subroutine finalize
	use etc,            only: istep,istart
	use mpi_params,     only: myid
	implicit none
	if(myid==0) &
		write(0,*) 'smooth exit after ',istep-istart-1,' time steps'
		call quit_mpi
		stop
end subroutine finalize


subroutine dinit(n,alpha,y)
	implicit none
	integer       :: n
	real(kind=8)  :: y(n),alpha
	y(:)=alpha
end subroutine dinit


subroutine LogMessage(msg,logfile)
	character(len=80)  :: msg,logfile
	open(1,file=logfile,position='append')
	write(1,*) msg
	close(1)
end subroutine LogMessage

real(kind=8) function myexp(arg)
	implicit none
	real(kind=8)    :: argmin = -64   ! exp(-256) = 6.616E-112  exp(-64)=1.603810890548638e-28 good enough
	real(kind=8)    :: argmax =  64   ! exp(256)  = 1.511E+111
	real(kind=8)    :: arg  ! assume negative value
	
	if( arg > 0 ) stop 'myexp expecting negative argument '
	if( arg < argmin ) then
		myexp = 0.d0
		return
	else
		myexp = exp(arg)
		return
	endif
end function myexp

real(kind=8) function myerf(arg)
	implicit none
	real(kind=8)    :: arg
	real(kind=8)    :: argmin = -10   ! ERF(arg) = -1.d0 for arg < -8, crash when |arg|>27
	real(kind=8)    :: argmax =  10   ! ERF(arg) =  1.d0 for arg >  8
	
	if( arg < argmin ) then
		myerf = -1.d0
		return
	elseif( arg > argmax) then
		myerf = 1.d0
		return
	else
		myerf = ERF(arg)
	endif
end function myerf

real(kind=8) function mytanh(arg)
	implicit none
	real(kind=8)    :: arg
	real(kind=8)    :: argmin = -10   
	real(kind=8)    :: argmax =  10   
	
	if( arg < argmin ) then
		mytanh = -1.d0
		return
	elseif( arg > argmax) then
		mytanh = 1.d0
		return
	else
		mytanh = tanh(arg)
	endif
end function mytanh



		
