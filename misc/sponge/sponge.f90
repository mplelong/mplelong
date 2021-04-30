!!=============================================
!!   Module for sponge regions to "soak" up
!!   motions at scales smaller than resolved in
!!   the parent simulation.
!!   to access variables, routines etc
!!   ==>   use sponge
!!=============================================
module sponge

	!-------------------------------------------------------------------------------------
	! Parameters defining the parent solution resolution
	!-------------------------------------------------------------------------------------
	real(kind=8)                     :: dx_parent = 4*(30000./512)     !! [m]
	real(kind=8)                     :: dy_parent = 1.d0           !! [m]
	real(kind=8)                     :: dz_parent = 4*(3000./256)    !! [m]
	real(kind=8), allocatable        :: sponge_filter_x(:)         !! [1] low pass x wavenumber filter
	real(kind=8), allocatable        :: sponge_filter_y(:)         !! [1] low pass y wavenumber filter
	real(kind=8), allocatable        :: sponge_filter_z(:)         !! [1] low pass z wavenumber filter
	real(kind=8), allocatable        :: FW(:,:,:)                  !! [1] filtered field weight (y,x,z)
	real(kind=8)                     :: offset = 0.025d0           !! fraction of domain to offset sponge region from bdries
	real(kind=8)                     :: width = 0.05d0             !! fraction of domain for sponge region width
									 !! offset and width generally reset in user_params module
	

contains

	subroutine setup_sponge_layers
		use mpi_params,               only: myid
		use independent_variables,    only: nx,ny,nz,Lx,Ly,Lz,x,y,z
		use differentiation_params,   only: kx,ky,kz
		use decomposition_params
		use methods_params,           only: do_sponging
		implicit none
		integer                          :: i,j,k,ig,kg
		real(kind=8)                     :: pi,beta,k_cut,dk
		real(kind=8)                     :: a,b,c
		real(kind=8), external           :: myexp,mytanh
		
		if( .NOT. do_sponging ) return
		
		allocate( sponge_filter_x(nx),sponge_filter_y(ny),sponge_filter_z(nz) )
		sponge_filter_x = 1.d0 ; sponge_filter_y = 1.d0 ; sponge_filter_z = 1.d0 
		
		
		pi = 4.d0*atan(1.d0)
		
		!----------------------
		! x wavenumber filter
		!----------------------
		if( nx > 1 ) then
			dk = pi/Lx
			beta = 2*dk
			k_cut = pi/dx_parent
		
			if( myid==0 ) open(1,file='output/debug_data/sponge_filter_x')
			do i=1,nx
				sponge_filter_x(i) = 0.5d0*( 1.d0 - mytanh((kx(i)-k_cut)/beta) )
				if(myid==0) write(1,*) x(i),sponge_filter_x(i)
			enddo
			if( myid==0 ) close(1)
		endif
		
		!----------------------
		! y wavenumber filter
		!----------------------
		if( ny > 1 ) then
			dk = pi/Ly
			beta = 2*dk
			k_cut = pi/dy_parent
		
			if( myid==0 ) open(1,file='output/debug_data/sponge_filter_y')
			do j=1,ny
				sponge_filter_y(j) = 0.5d0*( 1.d0 - mytanh((ky(j)-k_cut)/beta) )
				if(myid==0) write(1,*) y(j),sponge_filter_y(j)
			enddo
			if( myid==0 ) close(1)
		endif
		
		!----------------------
		! z wavenumber filter
		!----------------------
		if( nz > 1 ) then
			dk = pi/Lz
			beta = 2*dk
			k_cut = pi/dz_parent
		
			if( myid==0 ) open(1,file='output/debug_data/sponge_filter_z')
			do k=1,nz
				sponge_filter_z(k) = 0.5d0*( 1.d0 - mytanh((kz(k)-k_cut)/beta) )
				if(myid==0) write(1,*) z(k),sponge_filter_z(k)
			enddo
			if( myid==0 ) close(1)
		endif
		
		!---------------------------------------------------------------------
		! sponge region geometry, i.e. filter weighting as function of y,x,z
		!---------------------------------------------------------------------
		allocate( FW(ny,array_size(JDIM,YBLOCK,myid),array_size(KDIM,YBLOCK,myid)) )
		FW = 0.d0 
		do k=1,array_size(KDIM,YBLOCK,myid)
			kg = global_z_indices(START,YBLOCK,myid) + k - 1
			do i=1,array_size(JDIM,YBLOCK,myid)
				ig = global_x_indices(START,YBLOCK,myid) + i - 1
				do j=1,ny
				
				
					a = 0.d0  !sponge_weight(x(ig),Lx,offset,width)  ! 1d, depends on x
					b = 0.d0  !sponge_weight(y(j) ,Ly,offset,width)  ! 1d, depends on y
					c = sponge_weight(z(kg),Lz,offset,width)  ! 1d, depends on z
				
					!  close to bottom and top
					if( nz > 1 ) then
						if( z(kg) < offset*Lz .OR. z(kg) > (1.d0-offset)*Lz ) then
							!a=0.d0
							!b=0.d0
						endif
					endif
				
					! close to S/N
					if( ny > 1 ) then
						if( y(j) < offset*Ly .OR. y(j) > (1.d0-offset)*Ly ) then
							!a=0.d0
							!c=0.d0
						endif
					endif
					
					! close to E/W
					if( nx > 1 ) then
						if( x(ig) < offset*Lx .OR. x(ig) > (1.d0-offset)*Lx ) then
							!b=0.d0
							!c=0.d0
						endif
					endif
				
					! FW = 1 in region slightly inset from boundaries --> 0 at bdries and in interior
					FW(j,i,k) = MAXVAL((/a,b,c/))
													
				
				enddo
			enddo
		enddo
		if(myid==0) write(0,*) '................   Finished setting up sponge layers'
	 return
	end subroutine setup_sponge_layers
	
	real(kind=8) function tophat(x,x0,x1,beta_0,beta_1)
	!--------------------------------------------------------
	! smooth tophat function with different decay rates 
	! at left and right ends
	!--------------------------------------------------------
		implicit none
		real(kind=8)            :: x,x0,x1,beta_0,beta_1,ans
		real(kind=8), external  :: myexp
		integer, parameter      :: p=2
		if( x >= x0 .and. x <= x1 ) then
			ans = 1.d0
		elseif( x < x0 ) then
			ans = myexp(-((x-x0)/beta_0)**p)
		elseif( x > x1 ) then
			ans = myexp(-((x-x1)/beta_1)**p)
		endif
		tophat = ans
	return 
	end function tophat
	
	real(kind=8) function sponge_weight(x,Lx,offset,width)
	!--------------------------------------------------------
	! evaluates the tophat function to define a sponge_weight
	! value based on distance from end points
	!--------------------------------------------------------
		implicit none
		real(kind=8)            :: x,Lx,offset,width
		real(kind=8)            :: x0,x1,beta_0,beta_1,a1,a2
		real(kind=8), external  :: myexp
		
		x0 = 0.  ; x1 = width*Lx
		beta_0 = width*Lx  ; beta_1 = width*Lx
		a1 = tophat(x,x0,x1,beta_0,beta_1)
		
		x1 = Lx ; x0 = Lx - width*Lx
		beta_1 = width*Lx ; beta_0 = width*Lx
		a2 =  tophat(x,x0,x1,beta_0,beta_1)
					
		sponge_weight = a1 + a2		
	end function sponge_weight
	
end module sponge
