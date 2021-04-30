subroutine transform_xy(in,out,dir,exp_type)
!------------------------------------------------------
! Transform ENTIRE 3d YBLOCK data in x and y
! dir indicates forward or inverse transforms
! in and out can be the same array 
!------------------------------------------------------
	use mpi_params,             only: myid
	use decomposition_params
	use differentiation_params
	use etc
	use independent_variables,  only: nx,ny 
	use intermediate_variables, only: tmpX,tmpY
 
	implicit none     
	real(kind=8)                  ::  in( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid) )
 
	real(kind=8)                  :: out( array_size(IDIM,YBLOCK,myid),   &
                                          array_size(JDIM,YBLOCK,myid),   &
                                          array_size(KDIM,YBLOCK,myid) )
  
	real(kind=8),save             :: xfac,yfac,normfac
	integer                       :: j,k,dim,dim_start
	integer                       :: dir,mydir,i0
	character(len=80)             :: exp_type(2)
 
 
	if( nx==1 ) then
		dim_start=2      ! don't need x transform plans if nx=1
	else
		dim_start=1
	endif
        
	do dim=dim_start,2  !! x & y dimensions	 
		! re-use the sin plan from deriv setup
		if( trim(exp_type(dim))=='sin' ) then
			xy_plan(dim,1)=sin_plan(dim)  ! for in dim   n-2 pt sin transform
			xy_plan(dim,2)=sin_plan(dim)  ! inv in dim   ""
  
		! re-use the cos plan from deriv setup
		elseif( trim(exp_type(dim))=='cos' ) then
			xy_plan(dim,1)=cos_plan(dim)  ! for in dim   n pt cos transform
			xy_plan(dim,2)=cos_plan(dim)  ! inv in dim   ""
  
		! re-use the fourier plan from deriv setup
		elseif( trim(exp_type(dim))=='fourier' ) then
			xy_plan(dim,1)=fourier_plan(dim,1)  ! for in dim
			xy_plan(dim,2)=fourier_plan(dim,2)  ! inv in dim
		endif  
	enddo
 
	!---------------------
	! x normalization
	!---------------------
	if( nx > 1 ) then 
		if( trim(exp_type(1)) == 'cos' ) then
			xfac = 1.d0/(2.d0*(nx-1.d0))
		elseif( trim(exp_type(1)) == 'sin' ) then
			xfac = 1.d0/(2.d0*(nx-1.d0))
		elseif( trim(exp_type(1)) == 'fourier' ) then
			xfac = 1.d0/dfloat(nx)
		endif   
	elseif(nx==1) then
		xfac=1.d0
	endif
 
	!---------------------
	! y normalization
	!---------------------
	if( ny > 1 ) then   
		if( trim(exp_type(2)) == 'cos' ) then
			yfac = 1.d0/(2.d0*(ny-1.d0))
		elseif( trim(exp_type(2)) == 'sin' ) then
			yfac = 1.d0/(2.d0*(ny-1.d0))
		elseif( trim(exp_type(2)) == 'fourier' ) then
			yfac = 1.d0/dfloat(ny)
		endif   
	elseif( ny == 1 ) then
		yfac=1.d0
	endif
	normfac = xfac*yfac
 
	!! array index corresponding to forward/inverse transform direction
	if(dir==1) then
		mydir=1
	elseif(dir==-1) then
		mydir=2
	else
		stop ' dir is wrong in xy transform'
	endif
 
 
	!-------------------------------------------------------
	! do the x transforms
	!-------------------------------------------------------
	if( trim(exp_type(1)) == 'sin' ) then
		i0=2    ! exclude data endpoints for sin transform
	else
		i0=1
	endif
 
	if( nx > 1 ) then   
		!----------------------------------------
		!  transpose to XBLOCK format
		!----------------------------------------
		call yblock_2_xblock(in,tmpX)
  
		!------------------------------
		!  loop and x transform
		!------------------------------
		do k=1,array_size(KDIM,XBLOCK,myid)  
			do j=1,array_size(JDIM,XBLOCK,myid)
				call dfftw_execute_r2r(xy_plan(1,mydir),tmpX(i0,j,k,1),tmpX(i0,j,k,2))
			enddo
		enddo
   
		if( trim(exp_type(1)) == 'sin' ) then
			tmpX(1,:,:,2)=0.d0    ! "mean"/endpoint location in extended array
			tmpX(nx,:,:,2)=0.d0   ! "nyquist"/endpoint location in extended array
		endif

		!----------------------------------------
		!  transpose back to YBLOCK format
		!  store result in tmpY(:,:,:,6)
		!----------------------------------------
		if( ny > 1 ) then      ! result goes in tmpY
			call xblock_2_yblock(tmpX(1,1,1,2),tmpY(1,1,1,6))
		elseif( ny == 1 ) then  ! result goes directly into 'out'
			call xblock_2_yblock(tmpX(1,1,1,2),out)
			goto 999    ! don't do y transforms
		endif
  
 	endif
 

	!-------------------------------------------
	! y transform the results in tmpY(:,:,:,6)
	! and store the results in out
	!-------------------------------------------
	if( trim(exp_type(2)) == 'sin'  ) then
		i0=2    ! exclude data endpoints for sin transform
	else
		i0=1
	endif
   
	if( nx > 1 ) then       ! input is in tmpY 

		do k=1,array_size(KDIM,YBLOCK,myid)  
			do j=1,array_size(JDIM,YBLOCK,myid)
				call dfftw_execute_r2r(xy_plan(2,mydir),tmpY(i0,j,k,6),out(i0,j,k) )                                              
			enddo
		enddo
      
	elseif( nx==1 ) then   ! input still in 'in'
		do k=1,array_size(KDIM,YBLOCK,myid)  
			do j=1,array_size(JDIM,YBLOCK,myid)
				call dfftw_execute_r2r(xy_plan(2,mydir),in(i0,j,k),out(i0,j,k) )
			enddo 
		enddo
      
	endif
  
	if( trim(exp_type(2)) == 'sin' ) then
		out(1,:,:)=0.d0    ! "mean"/endpoint location in extended array
		out(ny,:,:)=0.d0   ! "nyquist"/endpoint location in extended array
	endif

999 continue
	if( dir==-1 ) then
		!-------------------------------------------
		! only normalize the inverse transforms
		!-------------------------------------------
		out = out*normfac
	endif 

 return
end subroutine transform_xy




subroutine test_transform_xy
	!------------------------------------------------------
	! Transform ENTIRE 3d YBLOCK data in x and y
	! dir indicates forward or inverse transforms
	! in and out can be the same array 
	!------------------------------------------------------
	use mpi_params,             only: myid
	use decomposition_params
	use intermediate_variables, only: tmpY
	use dimensional_scales,     only: length_scale
	use independent_variables,  only: x,y,Lx,Ly
	use etc
 
	implicit none     
	integer                       :: dir,i,j,k
	character(len=80)             :: exp_type(2)
	real(kind=8)                  :: pi,kx,ky,ans,diff,tol=1.e-15
 
 	!---------------------------------------------------------------
	!   test the 2d transform of a known function of x and y
	!---------------------------------------------------------------
	pi=4.*datan(1.d0)
	kx=1.d0*pi/Lx 
	ky=1.d0*pi/Ly   
 
	do k=1,array_size(KDIM,YBLOCK,myid)
		do i=1,array_size(JDIM,YBLOCK,myid)
			do j=1,array_size(IDIM,YBLOCK,myid)  
				tmpY(j,i,k,1) = cos(kx*x(i)) + sin(ky*y(j))     ! in
			enddo
		enddo
	enddo
 
 
	dir = 1   ! 1=for  -1=inv  
	exp_type(1)='cos'   ! x dir
	exp_type(2)='cos'   ! y dir
	call transform_xy(tmpY(1,1,1,1),tmpY(1,1,1,2),dir,exp_type)   ! in,out,dir...
  
	dir=-1
	call transform_xy(tmpY(1,1,1,2),tmpY(1,1,1,1),dir,exp_type)   ! in,out,dir...
 
	k=(array_size(KDIM,YBLOCK,myid)-1)/2   ! pick a "random" z level to test
	i=1
	do k=1,array_size(KDIM,YBLOCK,myid)
		do j=1,array_size(IDIM,YBLOCK,myid)
			ans = (cos(kx*x(i)) + sin(ky*y(j)))
			diff = abs( tmpY(j,i,k,1) - ans )
			if( diff .gt. tol ) then
				write(0,*) 'problem in transform_xy recovering test function ',j,i,k,myid
				stop
			endif
		enddo
	enddo
	
 	if(myid==0) then
  		message = '...  Test of FOR/INV transform_xy routine successful ...  out=in tol=1.d-15'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif 
 
 return
end subroutine test_transform_xy


