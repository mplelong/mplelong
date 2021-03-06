!-------------------------------------------------------------------------
! Routines to deal with boundary conditions. All these routines operate
! on boundary data already stored in arrays like east_vals(y,z,id) and 
! east_derivs(y,z,id) etc. These routines need not know how the data got 
! to those arrays. Filling these arrays at each time step is done in
! user_bvals_EW user_bvals_NS and user_bvals_BT at the moment.
!-------------------------------------------------------------------------
subroutine apply_bcs
	use mpi_params,                           only: myid
	use dependent_variables,                  only: u,v,w,s1,s2
	use independent_variables,                only: Lx,Ly,Lz,tnp1
	use methods_params,                       only: do_second_scalar,user_bcs
	use etc,                                  only: istep,istart
	use decomposition_params 
	implicit none
	integer, save                                :: locnx, ny, locnz, npts(3)
	real(kind=8),allocatable,dimension(:),save   :: x, y, z   ! local portions
	integer                                      :: i, j, k
	character(len=80)                            :: dir='xyz'
	logical,save                                 :: first_entry=.TRUE.
 
 	if( .NOT. user_bcs ) return    ! some configurations never need bc info
 	
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
		
		!-----------------------------------
		! smoothing distance in grid points
		!-----------------------------------
		npts(:)=6
		
		first_entry=.FALSE.
	endif
	
	dir = 'xyz'
	call boundary_smooth(u,dir,npts)
	call boundary_smooth(v,dir,npts)
	call boundary_smooth(w,dir,npts)
	call boundary_smooth(s1,dir,npts)
	
	!--------------------------------------------------------
	! extrapolate boundary info about tangential velocity 
	! components across thin boundary layers where finite
	! width boundary jumps are applied (i.e. eqns are inexact)
	!--------------------------------------------------------
	if( mod(istep,2)==0 ) then
		call extrapolate_to_boundaries     ! info from interior to boundaries
	else
		call extrapolate_from_boundaries   ! info from bdries to interior	
	endif
	
	
	!--------------------------------------------------------
	! below here is just testing output 
	!--------------------------------------------------------
	if( istep > istart ) return
	
	if( istep==istart .and. tnp1 == 0.d0 ) then   ! initial, immediate test applying BCs to ICS
	
		if(myid==0) then
 			open(1,file='output/debug_data/check_NS_bdry_scheme')
 				do j=1,ny
 					write(1,*) y(j),u(j,1,1),v(j,1,1),w(j,1,1),s1(j,1,1)
 				enddo	
 			close(1)
 		endif
 		
 		if( z(1)==0.d0 ) then
 			open(1,file='output/debug_data/check_B_bdry_scheme')			
 				do k=1,locnz
 					write(1,*) z(k),u(1,1,k),v(1,1,k),w(1,1,k),s1(1,1,k)
 				enddo
 			close(1)
 		endif
 		
 		if( z(locnz)==Lz ) then	 		
 			open(1,file='output/debug_data/check_T_bdry_scheme')
 				do k=1,locnz
 					write(1,*) z(k),u(1,1,k),v(1,1,k),w(1,1,k),s1(1,1,k)
 				enddo
 			close(1)
 		endif 
 		
 	endif
 	
			
	if( istep==istart .and. tnp1 > 0.d0) then	! test at end of first step
	
		if(myid==0) then
			i=1 ; k=locnz/2
 			open(1,file='output/debug_data/final_u_of_y')
 				do j=1,ny
 					write(1,*) y(j),u(j,i,k),v(j,i,k),w(j,i,k),s1(j,i,k)
 				enddo	
 			close(1)
 		endif		
 		if( z(1)==0.d0 ) then
 			open(1,file='output/debug_data/final_u_of_z_B')
 				i = 1; j = ny/2 + 1
 				do k=1,locnz
 					write(1,*) z(k),u(j,i,k),v(j,i,k),w(j,i,k),s1(j,i,k)
 				enddo
 			close(1)
 		endif 		 		
 		if( z(locnz)==Lz ) then
 			i = 1; j = ny/2 + 1
 			open(1,file='output/debug_data/final_u_of_z_T')
 				do k=1,locnz
 					write(1,*) z(k),u(j,i,k),v(j,i,k),w(j,i,k),s1(j,i,k)
 				enddo
 			close(1)
 		endif 	
 			 		
	endif
			
 return
end subroutine apply_bcs

subroutine extrapolate_from_boundaries
!----------------------------------------------------------------------
! Routine to extrapolate boundary knowledge about tangential velocity
! components into the interior within a thin boundary layer.
! Also insert scalar boundary values at domain edges.
!----------------------------------------------------------------------
	use mpi_params,            only: myid
	use decomposition_params
	use boundary_data
	use dependent_variables,   only: u,v,w,s1,s2
	use independent_variables, only: Lx,Ly,Lz,x_periodic,y_periodic,z_periodic,z_FSRL
	use methods_params,        only: do_second_scalar
	
	implicit none
	integer                        :: i,j,k,id
	integer, save                  :: ii=4, jj=4, kk=4   ! matching locations  (4 better than 6)
	integer, save                  :: locnx,ny,locnz
	real(kind=8),allocatable,save  :: x(:),y(:),z(:) 
	real(kind=8),allocatable,save  :: f_of_x(:),f_of_y(:),f_of_z(:)
	real(kind=8)                   :: offset
	logical, save                  :: first_entry=.TRUE.     
	
	if( first_entry ) then
		locnx = array_size(JDIM,YBLOCK,myid)
		ny  =   array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		
		allocate( x(locnx), y(ny), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)
		
		allocate( f_of_x(locnx), f_of_y(ny), f_of_z(locnz) )
		f_of_x = 0.d0 ; f_of_y = 0.d0 ; f_of_z = 0.d0 
		
		if( locnx > 1 .and. locnx <= ii ) stop 'locnx too small for fix_boundary_layers'
		if( locnz > 1 .and. locnz <= kk ) stop 'locnz too small for fix_boundary_layers'		
		first_entry = .FALSE.
	endif
	
	
	! fix BLs near east and west boundaries x=0 and x=Lx
	if( locnx > 1 .and. .NOT. x_periodic  ) then
		if( x(1)==0.d0 ) then     ! east bdry
			do k=1,locnz
				do j=1,ny
				
					id = 1  ! u
					f_of_x(1:ii) = east_vals(j,k,id) + east_derivs(j,k,id)*(x(1:ii)-0.d0)  ! 1st order Taylor series enough
					offset = u(j,ii,k) - f_of_x(ii)
					u(j,1:ii,k) = f_of_x(1:ii) + offset
					
					id = 2  ! v
					f_of_x(1:ii) = east_vals(j,k,id) + east_derivs(j,k,id)*(x(1:ii)-0.d0)  ! 1st order Taylor series enough
					offset = v(j,ii,k) - f_of_x(ii)
					v(j,1:ii,k) = f_of_x(1:ii) + offset
					
					id = 3  ! w
					f_of_x(1:ii) = east_vals(j,k,id) + east_derivs(j,k,id)*(x(1:ii)-0.d0)  ! 1st order Taylor series enough
					offset = w(j,ii,k) - f_of_x(ii)
					w(j,1:ii,k) = f_of_x(1:ii) + offset
					
					id = 4  ! s1
					f_of_x(1:ii) = east_vals(j,k,id) + east_derivs(j,k,id)*(x(1:ii)-0.d0)  ! 1st order Taylor series enough
					offset = s1(j,ii,k) - f_of_x(ii)
					s1(j,1:ii,k) = f_of_x(1:ii) + offset
					
					if( do_second_scalar) s2(j,1,k) = east_vals(j,k,5)
				enddo
			enddo			
		endif
		if( x(locnx)==Lx ) then   ! west bdry
			do k=1,locnz
				do j=1,ny
					id = 1  ! u
					f_of_x(locnx-ii+1:locnx) = west_vals(j,k,id) + west_derivs(j,k,id)*(x(locnx-ii+1:locnx)-Lx) 
					offset = u(j,locnx-ii+1,k) - f_of_x(locnx-ii+1)
					u(j,locnx-ii+1:locnx,k) = f_of_x(locnx-ii+1:locnx) + offset
					
					id = 2  ! v
					f_of_x(locnx-ii+1:locnx) = west_vals(j,k,id) + west_derivs(j,k,id)*(x(locnx-ii+1:locnx)-Lx) 
					offset = v(j,locnx-ii+1,k) - f_of_x(locnx-ii+1)
					v(j,locnx-ii+1:locnx,k) = f_of_x(locnx-ii+1:locnx) + offset
					
					id = 3  ! w
					f_of_x(locnx-ii+1:locnx) = west_vals(j,k,id) + west_derivs(j,k,id)*(x(locnx-ii+1:locnx)-Lx) 
					offset = w(j,locnx-ii+1,k) - f_of_x(locnx-ii+1)
					w(j,locnx-ii+1:locnx,k) = f_of_x(locnx-ii+1:locnx) + offset
					
					id = 4  ! s1
					f_of_x(locnx-ii+1:locnx) = west_vals(j,k,id) + west_derivs(j,k,id)*(x(locnx-ii+1:locnx)-Lx) 
					offset = s1(j,locnx-ii+1,k) - f_of_x(locnx-ii+1)
					s1(j,locnx-ii+1:locnx,k) = f_of_x(locnx-ii+1:locnx) + offset
					
					if( do_second_scalar) s2(j,locnx,k) = west_vals(j,k,5)
				enddo
			enddo
		endif
	endif
	
	! fix BLs near south and north boundaries y=0 and y=Ly
	if( ny > 1  .and. .NOT. y_periodic ) then
		do k=1,locnz
			do i=1,locnx
				id = 1  ! u
				f_of_y(1:jj) = south_vals(i,k,id) + south_derivs(i,k,id)*(y(1:jj)-0.d0)  ! 1st order Taylor series enough
				offset = u(jj,i,k) - f_of_y(jj)
				u(1:jj,i,k) = f_of_y(1:jj) + offset
				
				id = 2  ! v
				f_of_y(1:jj) = south_vals(i,k,id) + south_derivs(i,k,id)*(y(1:jj)-0.d0)  ! 1st order Taylor series enough
				offset = v(jj,i,k) - f_of_y(jj)
				v(1:jj,i,k) = f_of_y(1:jj) + offset
				
				id = 3  ! w
				f_of_y(1:jj) = south_vals(i,k,id) + south_derivs(i,k,id)*(y(1:jj)-0.d0)  ! 1st order Taylor series enough
				offset = w(jj,i,k) - f_of_y(jj)
				w(1:jj,i,k) = f_of_y(1:jj) + offset
				
				id = 4  ! s1
				f_of_y(1:jj) = south_vals(i,k,id) + south_derivs(i,k,id)*(y(1:jj)-0.d0)  ! 1st order Taylor series enough
				offset = s1(jj,i,k) - f_of_y(jj)
				s1(1:jj,i,k) = f_of_y(1:jj) + offset
				
				id = 1  ! u
				f_of_y(ny-jj+1:ny) = north_vals(i,k,id) + north_derivs(i,k,id)*(y(ny-jj+1:ny)-Ly)  ! 1st order Taylor series enough
				offset = u(ny-jj+1,i,k) - f_of_y(ny-jj+1)
				u(ny-jj+1:ny,i,k) = f_of_y(ny-jj+1:ny) + offset
				
				id = 2  ! v
				f_of_y(ny-jj+1:ny) = north_vals(i,k,id) + north_derivs(i,k,id)*(y(ny-jj+1:ny)-Ly)  ! 1st order Taylor series enough
				offset = v(ny-jj+1,i,k) - f_of_y(ny-jj+1)
				v(ny-jj+1:ny,i,k) = f_of_y(ny-jj+1:ny) + offset
				
				id = 3  ! w
				f_of_y(ny-jj+1:ny) = north_vals(i,k,id) + north_derivs(i,k,id)*(y(ny-jj+1:ny)-Ly)  ! 1st order Taylor series enough
				offset = w(ny-jj+1,i,k) - f_of_y(ny-jj+1)
				w(ny-jj+1:ny,i,k) = f_of_y(ny-jj+1:ny) + offset
				
				id = 4  ! s1
				f_of_y(ny-jj+1:ny) = north_vals(i,k,id) + north_derivs(i,k,id)*(y(ny-jj+1:ny)-Ly)  ! 1st order Taylor series enough
				offset = s1(ny-jj+1,i,k) - f_of_y(ny-jj+1)
				s1(ny-jj+1:ny,i,k) = f_of_y(ny-jj+1:ny) + offset
				
				if( do_second_scalar) s2(1,i,k) = south_vals(i,k,5)
				
				if( do_second_scalar) s2(ny,i,k) = north_vals(i,k,5)							
			enddo
		enddo
		
		
	endif
	
	! fix BLs near bottom and top boundaries
	if( locnz > 1  .and. .NOT. z_periodic .and. .NOT. z_FSRL ) then
		if( z(1)==0.d0 ) then     ! bottom bdry
			do i=1,locnx
				do j=1,ny
					id = 1  ! u
					f_of_z(1:kk) = bottom_vals(i,j,id) + bottom_derivs(i,j,id)*(z(1:kk)-0.d0)  ! 1st order Taylor series enough
					offset = u(j,i,kk) - f_of_z(kk)
					u(j,i,1:kk) = f_of_z(1:kk) + offset
					
					id = 2  ! v
					f_of_z(1:kk) = bottom_vals(i,j,id) + bottom_derivs(i,j,id)*(z(1:kk)-0.d0)  ! 1st order Taylor series enough
					offset = v(j,i,kk) - f_of_z(kk)
					v(j,i,1:kk) = f_of_z(1:kk) + offset
					
					id = 3  ! w
					f_of_z(1:kk) = bottom_vals(i,j,id) + bottom_derivs(i,j,id)*(z(1:kk)-0.d0)  ! 1st order Taylor series enough
					offset = w(j,i,kk) - f_of_z(kk)
					w(j,i,1:kk) = f_of_z(1:kk) + offset
					
					id = 4  ! s1
					f_of_z(1:kk) = bottom_vals(i,j,id) + bottom_derivs(i,j,id)*(z(1:kk)-0.d0)  ! 1st order Taylor series enough
					offset = s1(j,i,kk) - f_of_z(kk)
					s1(j,i,1:kk) = f_of_z(1:kk) + offset
					
					if( do_second_scalar) s2(j,i,1) = bottom_vals(i,j,5)
				enddo
			enddo			
		endif
		
		if( z(locnz)==Lz ) then   ! top bdry
			do i=1,locnx
				do j=1,ny
					id = 1  ! u
					f_of_z(locnz-kk+1:locnz) = top_vals(i,j,id) + top_derivs(i,j,id)*(z(locnz-kk+1:locnz)-Lz)  ! 1st order Taylor series enough
					offset = u(j,i,locnz-kk+1) - f_of_z(locnz-kk+1)
					u(j,i,locnz-kk+1:locnz) = f_of_z(locnz-kk+1:locnz) + offset
					
					id = 2  ! v
					f_of_z(locnz-kk+1:locnz) = top_vals(i,j,id) + top_derivs(i,j,id)*(z(locnz-kk+1:locnz)-Lz)  ! 1st order Taylor series enough
					offset = v(j,i,locnz-kk+1) - f_of_z(locnz-kk+1)
					v(j,i,locnz-kk+1:locnz) = f_of_z(locnz-kk+1:locnz) + offset
					
					id = 3  ! w
					f_of_z(locnz-kk+1:locnz) = top_vals(i,j,id) + top_derivs(i,j,id)*(z(locnz-kk+1:locnz)-Lz)  ! 1st order Taylor series enough
					offset = w(j,i,locnz-kk+1) - f_of_z(locnz-kk+1)
					w(j,i,locnz-kk+1:locnz) = f_of_z(locnz-kk+1:locnz) + offset
					
					id = 4  ! s1
					f_of_z(locnz-kk+1:locnz) = top_vals(i,j,id) + top_derivs(i,j,id)*(z(locnz-kk+1:locnz)-Lz)  ! 1st order Taylor series enough
					offset = s1(j,i,locnz-kk+1) - f_of_z(locnz-kk+1)
					s1(j,i,locnz-kk+1:locnz) = f_of_z(locnz-kk+1:locnz) + offset
										
					if( do_second_scalar) s2(j,i,locnz) = top_vals(i,j,5)
				enddo
			enddo
		endif
	endif
 return	
end subroutine extrapolate_from_boundaries


subroutine extrapolate_to_boundaries
!----------------------------------------------------------------------
! Routine to extrapolate solutions from just inside a boundary layer
! to the boundary.
! Also insert scalar boundary values at domain edges.
!----------------------------------------------------------------------
	use mpi_params,            only: myid
	use decomposition_params
	use boundary_data
	use dependent_variables,   only: u,v,w,s1,s2
	use independent_variables, only: Lx,Ly,Lz,x_periodic,y_periodic,z_periodic,z_FSRL
	use methods_params,        only: do_second_scalar
	
	implicit none
	integer                        :: i,j,k,id,iloc,jloc,kloc
	integer, save                  :: ii=6, jj=6, kk=6   ! expansion points
	integer, save                  :: locnx,ny,locnz
	real(kind=8),allocatable,save  :: x(:),y(:),z(:) 
	real(kind=8),allocatable,save  :: f_of_x(:),f_of_y(:),f_of_z(:)
	real(kind=8), save             :: dx,dy,dz
	real(kind=8)                   :: offset,m,x0,y0
	logical, save                  :: first_entry=.TRUE.     
	
	if( first_entry ) then
		locnx = array_size(JDIM,YBLOCK,myid)
		ny  =   array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		
		
		allocate( x(locnx), y(ny), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)
		
		if(locnx>1) dx = x(2)-x(1)
		if(ny>1)    dy = y(2)-y(1)
		if(locnz>1) dz = z(2)-z(1)
		
		allocate( f_of_x(locnx), f_of_y(ny), f_of_z(locnz) )
		f_of_x = 0.d0 ; f_of_y = 0.d0 ; f_of_z = 0.d0 
		
		if( locnx > 1 .and. locnx <= ii ) stop 'locnx too small for fix_boundary_layers'
		if( locnz > 1 .and. locnz <= kk ) stop 'locnz too small for fix_boundary_layers'		
		first_entry = .FALSE.
	endif
	
	
	! fix BLs near east and west boundaries x=0 and x=Lx
	if( locnx > 1 .and. .NOT. x_periodic ) then
		if( x(1)==0.d0 ) then     ! east bdry
			do k=1,locnz
				do j=1,ny
				
					!   adjust u near x=0
					iloc = ii
					x0 = x(iloc)
					
					m = (u(j,iloc,k)-u(j,iloc-1,k))/dx        ! slope
					y0 = u(j,iloc,k)
					u(j,1:iloc,k) = y0 + m*( x(1:iloc)-x0 )   ! pt slope formula
					
					!   adjust v near x=0
					m = (v(j,iloc,k)-v(j,iloc-1,k))/dx        ! slope
					y0 = v(j,iloc,k)
					v(j,1:iloc,k) = y0 + m*( x(1:iloc)-x0 )   ! pt slope formula
					
					!   adjust w near x=0
					m = (w(j,iloc,k)-w(j,iloc-1,k))/dx        ! slope
					y0 = w(j,iloc,k)
					w(j,1:iloc,k) = y0 + m*( x(1:iloc)-x0 )   ! pt slope formula
					
					!   adjust s1 near x=0
					m = (s1(j,iloc,k)-s1(j,iloc-1,k))/dx        ! slope
					y0 = s1(j,iloc,k)
					s1(j,1:iloc,k) = y0 + m*( x(1:iloc)-x0 )   ! pt slope formula
										
					!s1(j,1,k) = east_vals(j,k,4)
					if( do_second_scalar) s2(j,1,k) = east_vals(j,k,5)
					
				enddo
			enddo			
		endif
		if( x(locnx)==Lx ) then   ! west bdry
			do k=1,locnz
				do j=1,ny
				
					!   adjust u near x=Lx
					iloc = locnx - ii + 1
					x0 = x(iloc)
					
					m = (u(j,iloc,k)-u(j,iloc-1,k))/dx        ! slope
					y0 = u(j,iloc,k)
					u(j,iloc:locnx,k) = y0 + m*( x(iloc:locnx)-x0 ) ! pt slope formula	
					
					!   adjust v near x=Lx
					m = (v(j,iloc,k)-v(j,iloc-1,k))/dx        ! slope
					y0 = v(j,iloc,k)
					v(j,iloc:locnx,k) = y0 + m*( x(iloc:locnx)-x0 ) ! pt slope formula

					!   adjust w near x=Lx
					m = (w(j,iloc,k)-w(j,iloc-1,k))/dx        ! slope
					y0 = w(j,iloc,k)
					w(j,iloc:locnx,k) = y0 + m*( x(iloc:locnx)-x0 ) ! pt slope formula
					
					!   adjust s1 near x=Lx
					m = (s1(j,iloc,k)-s1(j,iloc-1,k))/dx        ! slope
					y0 = s1(j,iloc,k)
					s1(j,iloc:locnx,k) = y0 + m*( x(iloc:locnx)-x0 ) ! pt slope formula
															
					if( do_second_scalar) s2(j,locnx,k) = west_vals(j,k,5)
				enddo
			enddo
		endif
	endif
	
	! fix BLs near south and north boundaries y=0 and y=Ly
	if( ny > 1  .and. .NOT. y_periodic ) then
		do k=1,locnz
			do i=1,locnx
				
				! adjust u near y=0
				jloc = jj
				x0 = y(jloc)
				
				m = (u(jloc,i,k)-u(jloc-1,i,k))/dy        ! slope
				y0 = u(jloc,i,k)
				u(1:jloc,i,k) = y0 + m*( y(1:jloc)-x0 )   ! pt slope formula
				
				! adjust v near y=0	
				m = (v(jloc,i,k)-v(jloc-1,i,k))/dy        ! slope
				y0 = v(jloc,i,k)
				v(1:jloc,i,k) = y0 + m*( y(1:jloc)-x0 )   ! pt slope formula
				
				! adjust w near y=0	
				m = (w(jloc,i,k)-w(jloc-1,i,k))/dy        ! slope
				y0 = w(jloc,i,k)
				w(1:jloc,i,k) = y0 + m*( y(1:jloc)-x0 )   ! pt slope formula
				
				! adjust s1 near y=0	
				m = (s1(jloc,i,k)-s1(jloc-1,i,k))/dy      ! slope
				y0 = s1(jloc,i,k)
				s1(1:jloc,i,k) = y0 + m*( y(1:jloc)-x0 )  ! pt slope formula
				
				
				! adjust u near y=Ly
				jloc = ny - jj + 1
				x0 = y(jloc)
				
				m = (u(jloc,i,k)-u(jloc-1,i,k))/dy        ! slope
				y0 = u(jloc,i,k)
				u(jloc:ny,i,k) = y0 + m*( y(jloc:ny)-x0 ) ! pt slope formula
				
				! adjust v near y=0	
				m = (v(jloc,i,k)-v(jloc-1,i,k))/dy        ! slope
				y0 = v(jloc,i,k)
				v(jloc:ny,i,k) = y0 + m*( y(jloc:ny)-x0 )
				
				! adjust w near y=0	
				m = (w(jloc,i,k)-w(jloc-1,i,k))/dy        ! slope
				y0 = w(jloc,i,k)
				w(jloc:ny,i,k) = y0 + m*( y(jloc:ny)-x0 )
				
				! adjust s1 near y=0	
				m = (s1(jloc,i,k)-s1(jloc-1,i,k))/dy      ! slope
				y0 = s1(jloc,i,k)
				s1(jloc:ny,i,k) = y0 + m*( y(jloc:ny)-x0 )
				
				
				if( do_second_scalar) s2(1,i,k) = south_vals(i,k,5)
				
				if( do_second_scalar) s2(ny,i,k) = north_vals(i,k,5)							
			enddo
		enddo
		
		
	endif
	
	! fix BLs near bottom and top boundaries
	if( locnz > 1  .and. .NOT. z_periodic .and. .NOT. z_FSRL ) then
		if( z(1)==0.d0 ) then     ! bottom bdry
			do i=1,locnx
				do j=1,ny
				
					! adjust u near z=0
					kloc = kk
					x0 = z(kloc)
					
					m = (u(j,i,kloc)-u(j,i,kloc-1))/dz        ! slope
					y0 = u(j,i,kloc)
					u(j,i,1:kloc) = y0 + m*( z(1:kloc)-x0 )   ! pt slope formula
					
					! adjust v near z=0
					m = (v(j,i,kloc)-v(j,i,kloc-1))/dz        ! slope
					y0 = v(j,i,kloc)
					v(j,i,1:kloc) = y0 + m*( z(1:kloc)-x0 )   ! pt slope formula
					
					! adjust w near z=0
					m = (w(j,i,kloc)-w(j,i,kloc-1))/dz        ! slope
					y0 = w(j,i,kloc)
					w(j,i,1:kloc) = y0 + m*( z(1:kloc)-x0 )   ! pt slope formula
					
					! adjust s1 near z=0
					m = (s1(j,i,kloc)-s1(j,i,kloc-1))/dz        ! slope
					y0 = s1(j,i,kloc)
					s1(j,i,1:kloc) = y0 + m*( z(1:kloc)-x0 )   ! pt slope formula
									
					if( do_second_scalar) s2(j,i,1) = bottom_vals(i,j,5)
				enddo
			enddo			
		endif
		
		if( z(locnz)==Lz ) then   ! top bdry
			do i=1,locnx
				do j=1,ny
					! adjust u near z=0
					kloc = locnz - kk + 1
					x0 = z(kloc)
					
					m = (u(j,i,kloc)-u(j,i,kloc-1))/dz        ! slope
					y0 = u(j,i,kloc)
					u(j,i,kloc:locnz) = y0 + m*( z(kloc:locnz)-x0 ) ! pt slope formula
					
					! adjust v near z=0
					m = (v(j,i,kloc)-v(j,i,kloc-1))/dz        ! slope
					y0 = v(j,i,kloc)
					v(j,i,kloc:locnz) = y0 + m*( z(kloc:locnz)-x0 )
					
					! adjust w near z=0
					m = (w(j,i,kloc)-w(j,i,kloc-1))/dz        ! slope
					y0 = w(j,i,kloc)
					w(j,i,kloc:locnz) = y0 + m*( z(kloc:locnz)-x0 )
					
					! adjust s1 near z=0
					m = (s1(j,i,kloc)-s1(j,i,kloc-1))/dz        ! slope
					y0 = s1(j,i,kloc)
					s1(j,i,kloc:locnz) = y0 + m*( z(kloc:locnz)-x0 )
					
					if( do_second_scalar) s2(j,i,locnz) = top_vals(i,j,5)
				enddo
			enddo
		endif
	endif
	
 return	
end subroutine extrapolate_to_boundaries



subroutine boundary_smooth(f,dir,npts)
!------------------------------------------------------------------------
!    smooth f near the boundaries
!    f is data array in YBLOCK storage, i.e. f(y,x,z)
!    dir is a key indicating in which directions to smooth
!    npts is the number of gridpoints in the decay scale for smoothing
!    x,y,z are the local coords for the calling processor
!------------------------------------------------------------------------
	use mpi_params,                  only: myid
	use decomposition_params
	use independent_variables,       only: x,y,z,x_periodic,y_periodic,z_periodic,z_FSRL
	use methods_params,              only: endpoint_smoothing
	implicit none
	integer, intent(in)                 :: npts(3)  ! x,y,z dirs respectively
	real(kind=8), intent(inout)         :: f(array_size(IDIM,YBLOCK,myid),     &
                                             array_size(JDIM,YBLOCK,myid),     &
                                             array_size(KDIM,YBLOCK,myid) )
    real(kind=8),allocatable,save       :: xvals(:),yvals(:),zvals(:)
	character(len=80)                   :: dir,idir
	real(kind=8)                        :: dx=1.d0,dy=1.d0,dz=1.d0
	real(kind=8), save                  :: gamma_x, gamma_y, gamma_z
	integer, save                       :: locnx, locnz, ny
	integer                             :: i,j,k
	logical, save                       :: first_entry=.TRUE.
	
	if( .NOT. endpoint_smoothing ) return
	
	if( first_entry ) then
		locnx = array_size(JDIM,YBLOCK,myid)              ! sizes, loop endpoints
		locnz = array_size(KDIM,YBLOCK,myid)
		ny = array_size(IDIM,YBLOCK,myid)
		if(locnx>1) dx = x(2)-x(1)  
		if(ny>1)    dy = y(2)-y(1)  
		if(locnz>1) dz = z(2)-z(1)  
		gamma_x = npts(1)*dx                    
		gamma_y = npts(2)*dy
		gamma_z = npts(3)*dz                   
		allocate( xvals(array_size(JDIM,YBLOCK,myid)) )   ! local coord. arrays
		allocate( yvals(array_size(IDIM,YBLOCK,myid)) )
		allocate( zvals(array_size(KDIM,YBLOCK,myid)) )
		call get_my_xvals( xvals, YBLOCK, myid )
		call get_my_yvals( yvals, YBLOCK, myid )
		call get_my_zvals( zvals, YBLOCK, myid )
		first_entry=.FALSE.
	endif
	
	if( locnx > 1 .and. .NOT. x_periodic ) then
		if( dir=='x' .or. dir=='xz' .or. dir=='xy' .or. dir=='xyz' ) then
			idir='x'    
			do k=1,locnz
				do j=1,ny
					call smooth_near_boundary(f(j,:,k),xvals,locnx,gamma_x,idir)
				enddo
			enddo
		endif
	endif
	
	if( ny > 1 .and. .NOT. y_periodic ) then
		if( dir=='y' .or. dir=='xy' .or. dir=='yz' .or. dir=='xyz' ) then
			idir='y'
			do k=1,locnz
				do i=1,locnx
					call smooth_near_boundary(f(:,i,k),yvals,ny,gamma_y,idir)
				enddo
			enddo
		endif
	endif
	
	if( locnz > 1 .and. .NOT. z_periodic .and. .NOT. z_FSRL) then
		if( dir=='z' .or. dir=='xz' .or. dir=='yz' .or. dir=='xyz' ) then
			idir='z'
			do i=1,locnx
				do j=1,ny
					call smooth_near_boundary(f(j,i,:),zvals,locnz,gamma_z,idir)
				enddo
			enddo
		endif
	endif
 return
end subroutine boundary_smooth


subroutine smooth_near_boundary(f,x,n,gamma,dir)
!--------------------------------------------------------------------------------
!   Given a 1d pencil of data and the corresponding coordinate values, both
!   generally local subsets of a full 1d slice, do weighted near-neighbor
!   smoothing with weights decreasing away from the boundary
!   N.B.  if the smoothing distance gamma is comparable to the span of x
!   there will be a decoupling across processors. The idea is that smoothing
!   is limited to the processors containing the boundary points
!--------------------------------------------------------------------------------
	use mpi_params,                  only: myid     
	use decomposition_params,        only: IDIM,JDIM,KDIM,YBLOCK,array_size
	use independent_variables,       only: Lx,Ly,Lz
	implicit none
	real(kind=8)                        :: f(n)         ! local data pencil
	real(kind=8)                        :: x(n)         ! local "x" pencil
	character(len=80)                   :: dir
	integer                             :: i,j,n,nmax,nsmooths
	real(kind=8), allocatable, save     :: fs(:)
	real(kind=8)                        :: L,w,gamma
	real(kind=8), external              :: myexp
	integer, save                       :: locnx,ny,locnz,p
	logical                             :: periodic
	logical, save                       :: first_entry=.TRUE.
	
	if( first_entry ) then
	
		p = 2    ! power of exponential decay window
		
		locnx = array_size(JDIM,YBLOCK,myid)
		   ny = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		
		nmax = maxval( (/locnx,ny,locnz/) )
		allocate( fs(nmax) )         ! for sure big enough for all calls
		fs = 0.d0
		
		first_entry=.FALSE.
	endif
			
	
	if( n==1 ) return   ! don't smooth collapsed dimensions
	
	if( dir=='x') then
		nsmooths=1
		L = Lx	
	elseif( dir=='y') then
		nsmooths=1
		L=Ly  
	elseif( dir=='z') then
		nsmooths=1
		L=Lz
	else
		stop ' problem w/ dir specification in smooth_near_boundary '
	endif
	
	do j=1,nsmooths
		! end values
		fs(1) = (f(1)+f(2))/2.d0  
		fs(n) = (f(n)+f(n-1))/2.d0
	
		! interior values
		do i=2,n-1
			! neighbor weights, 1 at boundary --> 0 in interior
			w = myexp(-(x(i)/gamma)**p) + myexp(-((L-x(i))/gamma)**p)
			fs(i) = ( w*f(i-1) + f(i) + w*f(i+1)) / (2.d0*w + 1.d0)    
		enddo
	
	
		if( x(1) == 0.d0 ) then
			fs(1) = fs(1)
		else
			w = myexp(-(x(i)/gamma)**p) 
			fs(1) = (f(1) + w*f(2)) / (w+1.d0)	
		endif
	
		if( x(n) == L ) then
			fs(n) = fs(n)
		else
			w = myexp(-((L-x(n))/gamma)**p)
			fs(n) = (f(n) + w*f(n-1)) / (w+1.d0)
		endif
	
		! overwrite f with the smoothed version
		f(1:n) = fs(1:n)
	enddo
	
 return
end subroutine smooth_near_boundary

subroutine fill_boundary_arrays
	use mpi_params,              only: myid
	use decomposition_params
	use boundary_data
	use independent_variables,   only: Lx,Ly,Lz,tnp1,x_periodic,y_periodic,z_periodic,dt
	use methods_params,          only: do_second_scalar,user_bcs
	implicit none
	integer                         :: id
	integer,save                    :: locnx,ny,locnz,nvars=4
	real(kind=8),allocatable,save   :: x(:),y(:),z(:)
	real(kind=8)                    :: XVAL,YVAL,ZVAL,tval
	character(len=80)               :: dir,side
	logical,save                    :: first_entry=.TRUE.
	
	if( .NOT. user_bcs ) return
	
	if( first_entry ) then
	
		if( do_second_scalar )  nvars = 5
		
		locnx = array_size(JDIM,YBLOCK,myid)
		ny    = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		allocate( x(locnx), y(ny), z(locnz) )
		call get_my_xvals(x,YBLOCK,myid)
		call get_my_yvals(y,YBLOCK,myid)
		call get_my_zvals(z,YBLOCK,myid)
		
		allocate( east_vals(ny,locnz,4),west_vals(ny,locnz,4) ) 
		allocate( north_vals(locnx,locnz,4),south_vals(locnx,locnz,4) )
		allocate( bottom_vals(locnx,ny,4),top_vals(locnx,ny,4) )
		east_vals=0.d0; west_vals=0.d0; north_vals=0.d0; 
		south_vals=0.d0; bottom_vals=0.d0; top_vals=0.d0
		
		allocate( east_derivs(ny,locnz,4),west_derivs(ny,locnz,4) ) 
		allocate( north_derivs(locnx,locnz,4),south_derivs(locnx,locnz,4) )
		allocate( bottom_derivs(locnx,ny,4),top_derivs(locnx,ny,4) )
		east_derivs=0.d0; west_derivs=0.d0; north_derivs=0.d0; 
		south_derivs=0.d0; bottom_derivs=0.d0; top_derivs=0.d0
				
		first_entry=.FALSE.
	endif
	
	tval = tnp1   ! during a time step, we need BVALSand normal DERIVS at the next time step
	
		
	if( boundary_data_source == 'datafiles' ) then
			
		do id=1,nvars
			
			if(ny>1 .and. .NOT. y_periodic) then
				call bvals_NS(x,y,z,tval,id,locnx,ny,locnz)
			endif
			
			if(locnz>1 .and. .NOT. z_periodic) then
				call bvals_BT(x,y,z,tval,id,locnx,ny,locnz)
			endif
			
			if(locnx>1 .and. .NOT. x_periodic) then
				call bvals_EW(x,y,z,tval,id,locnx,ny,locnz)   ! do EW last, it increments read counters
			endif
							
		enddo
				
			
	elseif( boundary_data_source == 'user_functions' ) then
			
		do id=1,nvars
	
			if(locnx>1 .and. .NOT. x_periodic) then
				side='E'
				call user_bvals_EW(x,y,z,tval,id,east_vals(1,1,id),east_derivs(1,1,id),locnx,ny,locnz,side,Lx)
				side='W'
				call user_bvals_EW(x,y,z,tval,id,west_vals(1,1,id),west_derivs(1,1,id),locnx,ny,locnz,side,Lx)
			endif
		
			if(ny>1 .and. .NOT. y_periodic) then
				side='S'
				call user_bvals_NS(x,y,z,tval,id,south_vals(1,1,id),south_derivs(1,1,id),locnx,ny,locnz,side,Ly)
				side='N'
				call user_bvals_NS(x,y,z,tval,id,north_vals(1,1,id),north_derivs(1,1,id),locnx,ny,locnz,side,Ly)
			endif
		
			if(locnz>1 .and. .NOT. z_periodic) then
				side='B'
				call user_bvals_BT(x,y,z,tval,id,bottom_vals(1,1,id),bottom_derivs(1,1,id),locnx,ny,locnz,side,Lz)
				side='T'
				call user_bvals_BT(x,y,z,tval,id,top_vals(1,1,id),top_derivs(1,1,id),locnx,ny,locnz,side,Lz)
			endif
		
		enddo
		
	endif	
			
 return
end subroutine fill_boundary_arrays

subroutine bvals_EW(x,y,z,t,id,nx,ny,nz)
	!--------------------------------------------------------------------------
	! flow_solve routine to read the boundary data netcdf files and
	! fill the arrays east_vals(y,z,id) and east_derivs(y,z,id)
	! and  west_vals(y,z,id) and west_derivs(y,z,id)
	! at the EAST and WEST boundaries x=0 and x=Lx at
	! the y and z values requested at time t. All values in dimensional units.
	!
	! x,y,z and nx,ny,nz are all local to the processor
	!--------------------------------------------------------------------------
	use mpi_params,             only: myid,comm,ierr
	use methods_params,         only: do_second_scalar
	use decomposition_params,   only: IDIM,JDIM,KDIM,array_size,YBLOCK, &
                                      proc_row,proc_col,p1,p2
    use independent_variables, only : Lx,t0,tf  ! initial and final integration times
    use boundary_data		
	implicit none
	integer, intent(in)                    :: id                  ! 1,2,3,4,5 for u,v,w,s1,s2
	integer, intent(in)                    :: nx,ny,nz            ! size of local portions of domain
	real(kind=8), intent(in)               :: t                   ! time in [s], global size Lx in [m]
	real(kind=8), intent(in)               :: x(nx),y(ny),z(nz)   ! local values of x, y and z
	real(kind=8), allocatable, save        :: VALS(:,:,:,:,:)     ! (y,z,t,id,iside)
	real(kind=8), allocatable, save        :: DERIVS(:,:,:,:,:)   ! (y,z,t,id,iside)	
	real(kind=8), allocatable, save        :: TIME(:)             ! time in boundary data files
	integer                                :: iside,i0,i1
	real(kind=8)                           :: t_left,t_right,xx
	logical, save                          :: do_east,do_west,do_read
	
	!------------------------------------------------------------
 	! filename stuff...
 	!------------------------------------------------------------
	integer                         ::  ncid,varid,i_proc,j_proc,dimid
	character(len=4)                ::  cid
	character(len=3), save          ::  iproc,jproc
	character(len=80)               ::  ncfile,val_name,deriv_name,char_side
	integer, save                   ::  start(3),count(3),nvars=4,read_time=.TRUE.
	logical, save                   ::  first_entry=.TRUE.
	integer                         ::  j,k,icount,remaining
	include 'netcdf.inc'
	
	if( first_entry ) then
	
		if( do_second_scalar ) nvars = 5
		
		t_end = 0.d0
		start_slice = 1
		
		allocate(   VALS(ny,nz,slices_per_chunk,nvars,2) )  ! E/W (y,z,t,id,side)
		allocate( DERIVS(ny,nz,slices_per_chunk,nvars,2) )
		VALS = 0.d0 ; DERIVS = 0.d0	
		allocate( TIME(slices_per_chunk) )
		TIME = 0.d0	
		
		write(unit=iproc,fmt='(I3.3)') proc_row(YBLOCK,myid) ! 3 digit char string for iproc
   		write(unit=jproc,fmt='(I3.3)') proc_col(YBLOCK,myid) ! 3 digit char string for jproc
   		
   		!-----------------------------------------------------------
		! determine if this processor needs to fill tha data values
		!-----------------------------------------------------------
		do_east=.FALSE. ;  do_west=.FALSE.
		if( x(1) == 0.d0 ) do_east = .TRUE.
		if( x(nx) == Lx )  do_west = .TRUE.
		
		if( myid == 0 ) then
			open(1,file='output/debug_data/BC_data_EW_specs')
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) ' initial and final integration times [s]   :',t0,tf
				write(1,*) ' processor myid=0  iproc,jproc             :',iproc,' , ',jproc
				write(1,*) ' chunk size ny, nz, slices_per_chunk,nvars :',ny,nz,slices_per_chunk,nvars
				write(1,*) ' for this processor do_east,do_west        :',do_east,do_west
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) ' '
			close(1)	
		endif
   	
		first_entry=.FALSE.
	endif
	   
	if( id == 1 ) then
		val_name = 'u' ; deriv_name = 'u_x'
	elseif( id == 2 ) then
		val_name = 'v' ; deriv_name = 'v_x'
	elseif( id == 3 ) then
		val_name = 'w' ; deriv_name = 'w_x'
	elseif( id == 4 ) then
		val_name = 's1' ; deriv_name = 's1_x'
	elseif( id == 5 ) then
		val_name = 's2' ; deriv_name = 's2_x'
	endif
			
	!-----------------------------------------------------------
	! determine if a new data chunk needs to be read and if so,
	! read the new chunk of boundary data
	!-----------------------------------------------------------
	if( do_east .or. do_west ) then
		
		if( t >= t_end ) then    ! a new chunk of boundary data is needed, read it
	
			do iside = 1,2    ! E to W
   		
   				if( iside==1 .and. do_east ) then
   					if(proc_row(YBLOCK,myid) .NE. 0) stop 'problem detecting do_east in bvals_EW' 
   					ncfile = 'BVALS/east_'//iproc//'-'//jproc//'.nc'  ! i.e. east_000-013.nc
   					do_read = .TRUE.   					
   				elseif( iside==2 .and. do_west ) then
   					if(proc_row(YBLOCK,myid) .NE. p1-1) stop 'problem detecting do_west in bvals_EW'
   					ncfile = 'BVALS/west_'//iproc//'-'//jproc//'.nc'  ! i.e. west_000-013.nc
   					do_read = .TRUE.
   				else
   					do_read = .FALSE.
   				endif
   			
   				if( do_read ) then
   					!--------------------------------------
   					!  open the file, check for error
   					!--------------------------------------
   					ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
   					if(ierr.ne.NF_NOERR) then
    					write(0,*) '... ERROR OPENING NETCDF FILE: bvals_EW ',trim(ncfile)
    					write(0,*) '... myid, ierr ',myid, ierr
    					stop
   					endif
   			
   					!-----------------------------------------
   					!  determine total number of time slices
   					!-----------------------------------------
   					ierr = NF_INQ_DIMID(ncid, "timedimension", dimid)
					ierr = NF_INQ_DIMLEN(ncid,dimid,num_time_slices)
					if(ierr.ne.NF_NOERR) then
    					write(0,*) '... ERROR GETTING TIMEDIMENSION LENGTH: bvals_EW ',trim(ncfile)
    					write(0,*) '... myid, ierr ',myid, ierr
    					stop
   					endif
					if( num_time_slices < slices_per_chunk ) then
						slices_per_chunk = num_time_slices  ! don't read more slices than there are in the file
					endif
					
					!---------------------------------------------
   					!  set start and count to read desired data
   					!---------------------------------------------
   					remaining = num_time_slices - start_slice + 1
   					icount = min(remaining,slices_per_chunk)
					start=(/1,1,start_slice/)
   					count=(/ny,nz,icount/)
   					
			
					!-----------------------------------------
   					!  get start and end times for this chunk
   					!  set the time increment
   					!-----------------------------------------
					ierr = NF_INQ_VARID(ncid,'time',varid)
					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start_slice,icount,TIME)
					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR GET_VARA_DOUBLE -> time: ',ierr,start_slice,icount    					
    					stop
   					endif
					t_start = TIME(1)
					dt_bvals = TIME(2)-TIME(1)
					
					if( t < TIME(1) .or. t > TIME(icount) ) then
						if(myid==0) then
							write(0,*) id,myid," time range of new chunk does not span current time t"
							write(0,*) id,myid,"       t,TIME(1),TIME(icount): ",t,TIME(1),TIME(icount)
							write(0,*) id,myid,"       start_slice,icount,num_time_slices: ",start_slice,icount,num_time_slices
						endif	
					endif
												   			
   					!------------------------------------------------
   					!  extract variable id for VAL , check for error
   					!------------------------------------------------
   					ierr=NF_INQ_VARID(ncid,val_name,varid)
   					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR INQ_VARID -> val, id: ',id,ierr
    					stop
   					endif

   					!------------------------------------------------
   					!  read the corresponding variable (dimensional)
   					!------------------------------------------------
   					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start,count,VALS(1,1,1,id,iside)) 
   					if (ierr.ne.NF_NOERR) then
   						write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> val, id: ',id,ierr
    					stop
   					endif
   				
   					!---------------------------------------------------
   					!  extract variable id for DERIV , check for error
   					!---------------------------------------------------
   					ierr=NF_INQ_VARID(ncid,deriv_name,varid)
   					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR INQ_VARID -> deriv, id: ',id,ierr
    					stop
   					endif

   					!------------------------------------------------
   					!  read the corresponding variable (dimensional)
   					!------------------------------------------------
   					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start,count,DERIVS(1,1,1,id,iside)) 
   					if (ierr.ne.NF_NOERR) then
   						write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> deriv, id: ',id,ierr
    					stop
   					endif
   				
   					if( myid == 0 .and. id<9 ) then
   						open(1,file='output/debug_data/BC_data_EW_specs',position='append')
   						write(1,*) ' finished reading file, iside        :',trim(ncfile),iside
   						write(1,*) ' id and variable names read          :',id,trim(val_name),' , ',trim(deriv_name)
   						write(1,*) ' total number of time slices in file :',num_time_slices
   						write(1,*) ' slice range for this chunk          :',start_slice,start_slice + icount - 1 
   						write(1,*) ' time range for this chunk           :',TIME(1),TIME(icount)
   						write(1,*) ' time value triggering read          :',t
   						write(1,*) '  '
   						write(1,*) '  '
   						close(1)
   					endif
   				
   				endif ! do_read block
   			
   			enddo   ! end loop over 2 endpoints
   			
   			! wait until all variables have read their data before incrementing
			if( id == nvars ) then   
				!--------------------------------------------------------------------------
				! increment indices used for deciding which slices are in the next chunk
				! after reading last variable, only do this for bvals_EW assuming this
				! is the LAST of the bvals routines to be called 
				!--------------------------------------------------------------------------
				t_end = TIME(icount)                    ! ending time of last chunk read
				start_slice = start_slice + icount - 1  ! starting slice for the next chunk to be read
			endif
   			   		   			
		endif   ! end block for reading a new chunk of boundary data
		
		!-----------------------------------------------------------------------
		! a valid boundary data chunk is now available (new or old)
		! interpolate between saved BC values to get values at the current time
		!-----------------------------------------------------------------------
		
   		i0 = floor( (t-t_start)/dt_bvals ) + 1        ! index of saved time just lt t
   		i1 = minval( (/i0 + 1,slices_per_chunk/) )    ! index of saved time just gt t
   		t_left = t_start + (i0-1)*dt_bvals            ! saved time in [s] just lt t
   		t_right = t_left + dt_bvals                   ! saved time in [s] just gt t
   		xx = (t - t_left)/dt_bvals                    ! linear interpolation weight
   		   
   		!---------------
   		! SANITY CHECK
   		!---------------
   		if( t < t_left .OR. t > t_right ) then
    		if(myid==0) then
     			write(0,*) 'misalignment in BC time values: ',t,t_start,t_end,dt_bvals
     			write(0,*) 't_left, t_right ',t_left, t_right
     			write(0,*) 'i0, i1, xx ',i0,i1,xx
     			STOP
    		endif
   		endif   
   		if( xx < 0.d0 .OR. xx > 1.d0 ) then
    		if(myid==0) then
     			write(0,*) 'misalignment in BC time values: ',t,t_start,t_end
     			write(0,*) 't, t_left, t_right ',t, t_left, t_right
     			write(0,*) 'xx, dt_bvals  ',xx, dt_bvals
     			STOP
    		endif
   		endif
   		   		   		 
   		!-------------------------------------------------
   		! do the interpolation and store the value
   		!-------------------------------------------------
   		if( do_east ) then
   			iside = 1
   			east_vals(:,:,id)   =   VALS(:,:,i0,id,iside) + xx*(   VALS(:,:,i1,id,iside) -   VALS(:,:,i0,id,iside) )  ! time interpolated values
   			east_derivs(:,:,id) = DERIVS(:,:,i0,id,iside) + xx*( DERIVS(:,:,i1,id,iside) - DERIVS(:,:,i0,id,iside) )  ! time interpolated derivs
   		endif
   		
   		if( do_west ) then
   			iside = 2
   			west_vals(:,:,id)   =   VALS(:,:,i0,id,iside) + xx*(   VALS(:,:,i1,id,iside) -   VALS(:,:,i0,id,iside) )  ! time interpolated values
   			west_derivs(:,:,id) = DERIVS(:,:,i0,id,iside) + xx*( DERIVS(:,:,i1,id,iside) - DERIVS(:,:,i0,id,iside) )  ! time interpolated derivs
   		endif
   		   		
	endif   ! end block for processors that own EW boundaries
	
	!-------------------------------------------------
   	! workload can be different, not all processors
   	! have boundaries and do the reads, sync back up
   	!-------------------------------------------------		
	call mpi_barrier(comm,ierr)  
	
 return	          
end subroutine bvals_EW

subroutine bvals_NS(x,y,z,t,id,nx,ny,nz)
	!--------------------------------------------------------------------------
	! flow_solve routine to read the boundary data netcdf files and
	! fill the arrays south_vals(x,z,id) and south_derivs(x,z,id)
	! and  north_vals(x,z,id) and north_derivs(x,z,id)
	! at the SOUTH and NORTH boundaries y=0 and y=Ly at
	! the x and z values requested at time t. All values in dimensional units.
	!
	! x,y,z and nx,ny,nz are all local to the processor
	!--------------------------------------------------------------------------
	use mpi_params,             only: myid,comm,ierr
	use methods_params,         only: do_second_scalar
	use decomposition_params,   only: IDIM,JDIM,KDIM,array_size,YBLOCK, &
                                      proc_row,proc_col,p1,p2
    use independent_variables, only : Ly,t0,tf  ! initial and final integration times
    use boundary_data		
	implicit none
	integer, intent(in)                    :: id                  ! 1,2,3,4,5 for u,v,w,s1,s2
	integer, intent(in)                    :: nx,ny,nz            ! size of local portions of domain
	real(kind=8), intent(in)               :: t                   ! time in [s], global size Lx in [m]
	real(kind=8), intent(in)               :: x(nx),y(ny),z(nz)   ! local values of x, y and z
	real(kind=8), allocatable, save        :: VALS(:,:,:,:,:)     ! (y,z,t,id,iside)
	real(kind=8), allocatable, save        :: DERIVS(:,:,:,:,:)   ! (y,z,t,id,iside)	
	real(kind=8), allocatable, save        :: TIME(:)             ! time in boundary data files
	integer                                :: iside,i0,i1
	real(kind=8)                           :: t_left,t_right,xx
	logical, save                          :: do_south,do_north,do_read
	
	!------------------------------------------------------------
 	! filename stuff...
 	!------------------------------------------------------------
	integer                         ::  ncid,varid,i_proc,j_proc,dimid
	character(len=4)                ::  cid
	character(len=3), save          ::  iproc,jproc
	character(len=80)               ::  ncfile,val_name,deriv_name,char_side
	integer, save                   ::  start(3),count(3),nvars=4,read_time=.TRUE.
	logical, save                   ::  first_entry=.TRUE.
	integer                         ::  i,k,icount,remaining
	include 'netcdf.inc'
	
	if( first_entry ) then
	
		if( do_second_scalar ) nvars = 5
		
		t_end = 0.d0
		start_slice = 1
		
		allocate(   VALS(nx,nz,slices_per_chunk,nvars,2) )  ! S/N (x,z,t,id,side)
		allocate( DERIVS(nx,nz,slices_per_chunk,nvars,2) )
		VALS = 0.d0 ; DERIVS = 0.d0	
		allocate( TIME(slices_per_chunk) )
		TIME = 0.d0	
		
		write(unit=iproc,fmt='(I3.3)') proc_row(YBLOCK,myid) ! 3 digit char string for iproc
   		write(unit=jproc,fmt='(I3.3)') proc_col(YBLOCK,myid) ! 3 digit char string for jproc
   		
   		!-----------------------------------------------------------
		! determine if this processor needs to fill tha data values
		!-----------------------------------------------------------
		do_south=.FALSE. ;  do_north=.FALSE.
		if( y(1) == 0.d0 ) do_south = .TRUE.
		if( y(ny) == Ly )  do_north = .TRUE.
		
		if( myid == 0 ) then
			open(1,file='output/debug_data/BC_data_NS_specs')
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) ' initial and final integration times [s]   :',t0,tf
				write(1,*) ' processor myid=0  iproc,jproc             :',iproc,' , ',jproc
				write(1,*) ' chunk size nx, nz, slices_per_chunk,nvars :',nx,nz,slices_per_chunk,nvars
				write(1,*) ' for this processor do_south,do_north      :',do_south,do_north
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) ' '
			close(1)	
		endif
   	
		first_entry=.FALSE.
	endif
	   
	if( id == 1 ) then
		val_name = 'u' ; deriv_name = 'u_y'
	elseif( id == 2 ) then
		val_name = 'v' ; deriv_name = 'v_y'
	elseif( id == 3 ) then
		val_name = 'w' ; deriv_name = 'w_y'
	elseif( id == 4 ) then
		val_name = 's1' ; deriv_name = 's1_y'
	elseif( id == 5 ) then
		val_name = 's2' ; deriv_name = 's2_y'
	endif
			
	!-----------------------------------------------------------
	! determine if a new data chunk needs to be read and if so,
	! read the new chunk of boundary data
	!-----------------------------------------------------------
	if( do_south .or. do_north ) then
		
		if( t >= t_end ) then    ! a new chunk of boundary data is needed, read it
	
			do iside = 1,2    ! S to N
   		
   				if( iside==1 .and. do_south ) then
   					ncfile = 'BVALS/south_'//iproc//'-'//jproc//'.nc'  ! i.e. south_000-013.nc
   					do_read = .TRUE.   					
   				elseif( iside==2 .and. do_north ) then
   					ncfile = 'BVALS/north_'//iproc//'-'//jproc//'.nc'  ! i.e. north_000-013.nc
   					do_read = .TRUE.
   				else
   					do_read = .FALSE.
   				endif
   			
   				if( do_read ) then
   					!--------------------------------------
   					!  open the file, check for error
   					!--------------------------------------
   					ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
   					if(ierr.ne.NF_NOERR) then
    					write(0,*) '... ERROR OPENING NETCDF FILE: bvals_NS ',trim(ncfile)
    					write(0,*) '... myid, ierr ',myid, ierr
    					stop
   					endif
   			
   					!-----------------------------------------
   					!  determine total number of time slices
   					!-----------------------------------------
   					ierr = NF_INQ_DIMID(ncid, "timedimension", dimid)
					ierr = NF_INQ_DIMLEN(ncid,dimid,num_time_slices)
					if(ierr.ne.NF_NOERR) then
    					write(0,*) '... ERROR GETTING TIMEDIMENSION LENGTH: bvals_NS ',trim(ncfile)
    					write(0,*) '... myid, ierr ',myid, ierr
    					stop
   					endif
					if( num_time_slices < slices_per_chunk ) then
						slices_per_chunk = num_time_slices  ! don't read more slices than there are in the file
					endif
					
					!---------------------------------------------
   					!  set start and count to read desired data
   					!---------------------------------------------
   					remaining = num_time_slices - start_slice + 1
   					icount = min(remaining,slices_per_chunk)
					start=(/1,1,start_slice/)
   					count=(/nx,nz,icount/)
   					
			
					!-----------------------------------------
   					!  get start and end times for this chunk
   					!  set the time increment
   					!-----------------------------------------
					ierr = NF_INQ_VARID(ncid,'time',varid)
					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start_slice,icount,TIME)
					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR GET_VARA_DOUBLE -> time: ',ierr,start_slice,icount    					
    					stop
   					endif
					t_start = TIME(1)
					dt_bvals = TIME(2)-TIME(1)
					
					if( t < TIME(1) .or. t > TIME(icount) ) then
						if(myid==0) then
							write(0,*) id,myid," time range of new chunk does not span current time t"
							write(0,*) id,myid,"       t,TIME(1),TIME(icount): ",t,TIME(1),TIME(icount)
							write(0,*) id,myid,"       start_slice,icount,num_time_slices: ",start_slice,icount,num_time_slices
						endif	
					endif
												   			
   					!------------------------------------------------
   					!  extract variable id for VAL , check for error
   					!------------------------------------------------
   					ierr=NF_INQ_VARID(ncid,val_name,varid)
   					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR INQ_VARID -> val, id: ',id,ierr
    					stop
   					endif

   					!------------------------------------------------
   					!  read the corresponding variable (dimensional)
   					!------------------------------------------------
   					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start,count,VALS(1,1,1,id,iside)) 
   					if (ierr.ne.NF_NOERR) then
   						write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> val, id: ',id,ierr
    					stop
   					endif
   				
   					!---------------------------------------------------
   					!  extract variable id for DERIV , check for error
   					!---------------------------------------------------
   					ierr=NF_INQ_VARID(ncid,deriv_name,varid)
   					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR INQ_VARID -> deriv, id: ',id,ierr
    					stop
   					endif

   					!------------------------------------------------
   					!  read the corresponding variable (dimensional)
   					!------------------------------------------------
   					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start,count,DERIVS(1,1,1,id,iside)) 
   					if (ierr.ne.NF_NOERR) then
   						write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> deriv, id: ',id,ierr
    					stop
   					endif
   				
   					if( myid == 0 .and. id<9 ) then
   						open(1,file='output/debug_data/BC_data_NS_specs',position='append')
   						write(1,*) ' finished reading file, iside        :',trim(ncfile),iside
   						write(1,*) ' id and variable names read          :',id,trim(val_name),' , ',trim(deriv_name)
   						write(1,*) ' total number of time slices in file :',num_time_slices
   						write(1,*) ' slice range for this chunk          :',start_slice,start_slice + icount - 1 
   						write(1,*) ' time range for this chunk           :',TIME(1),TIME(icount)
   						write(1,*) ' time value triggering read          :',t
   						write(1,*) '  '
   						write(1,*) '  '
   						close(1)
   					endif
   				
   				endif ! do_read block
   			
   			enddo   ! end loop over 2 endpoints
   			   			   		   			
		endif   ! end block for reading a new chunk of boundary data
		
		!-----------------------------------------------------------------------
		! a valid boundary data chunk is now available (new or old)
		! interpolate between saved BC values to get values at the current time
		!-----------------------------------------------------------------------
		
   		i0 = floor( (t-t_start)/dt_bvals ) + 1        ! index of saved time just lt t
   		i1 = minval( (/i0 + 1,slices_per_chunk/) )    ! index of saved time just gt t
   		t_left = t_start + (i0-1)*dt_bvals            ! saved time in [s] just lt t
   		t_right = t_left + dt_bvals                   ! saved time in [s] just gt t
   		xx = (t - t_left)/dt_bvals                    ! linear interpolation weight
   		   
   		!---------------
   		! SANITY CHECK
   		!---------------
   		if( t < t_left .OR. t > t_right ) then
    		if(myid==0) then
     			write(0,*) 'misalignment in BC time values: ',t,t_start,t_end,dt_bvals
     			write(0,*) 't_left, t_right ',t_left, t_right
     			write(0,*) 'i0, i1, xx ',i0,i1,xx
     			STOP
    		endif
   		endif   
   		if( xx < 0.d0 .OR. xx > 1.d0 ) then
    		if(myid==0) then
     			write(0,*) 'misalignment in BC time values: ',t,t_start,t_end
     			write(0,*) 't, t_left, t_right ',t, t_left, t_right
     			write(0,*) 'xx, dt_bvals  ',xx, dt_bvals
     			STOP
    		endif
   		endif
   		
   		!-------------------------------------------------
   		! do the interpolation and store the value
   		!-------------------------------------------------
   		if( do_south ) then
   			iside = 1
   			south_vals(:,:,id)   =   VALS(:,:,i0,id,iside) + xx*(   VALS(:,:,i1,id,iside) -   VALS(:,:,i0,id,iside) )  ! time interpolated values
   			south_derivs(:,:,id) = DERIVS(:,:,i0,id,iside) + xx*( DERIVS(:,:,i1,id,iside) - DERIVS(:,:,i0,id,iside) )  ! time interpolated derivs
   		endif
   		
   		if( do_north ) then
   			iside = 2
   			north_vals(:,:,id)   =   VALS(:,:,i0,id,iside) + xx*(   VALS(:,:,i1,id,iside) -   VALS(:,:,i0,id,iside) )  ! time interpolated values
   			north_derivs(:,:,id) = DERIVS(:,:,i0,id,iside) + xx*( DERIVS(:,:,i1,id,iside) - DERIVS(:,:,i0,id,iside) )  ! time interpolated derivs
   		endif
   		   		
	endif   ! end block for processors that own NS boundaries
	
	!-------------------------------------------------
   	! workload can be different, not all processors
   	! have boundaries and do the reads, sync back up
   	!-------------------------------------------------		
	call mpi_barrier(comm,ierr)  
	
 return	          
end subroutine bvals_NS

subroutine bvals_BT(x,y,z,t,id,nx,ny,nz)
	!--------------------------------------------------------------------------
	! flow_solve routine to read the boundary data netcdf files and
	! fill the arrays bottom_vals(x,y,id) and bottom_derivs(x,y,id)
	! and  top_vals(x,y,id) and top_derivs(x,y,id)
	! at the BOTTOM and TOP boundaries z=0 and z=Lz at
	! the x and y values requested at time t. All values in dimensional units.
	!
	! x,y,z and nx,ny,nz are all local to the processor
	!--------------------------------------------------------------------------
	use mpi_params,             only: myid,comm,ierr
	use methods_params,         only: do_second_scalar
	use decomposition_params,   only: IDIM,JDIM,KDIM,array_size,YBLOCK, &
                                      proc_row,proc_col,p1,p2
    use independent_variables, only : Lz,t0,tf  ! initial and final integration times
    use boundary_data		
	implicit none
	integer, intent(in)                    :: id                  ! 1,2,3,4,5 for u,v,w,s1,s2
	integer, intent(in)                    :: nx,ny,nz            ! size of local portions of domain
	real(kind=8), intent(in)               :: t                   ! time in [s], global size Lx in [m]
	real(kind=8), intent(in)               :: x(nx),y(ny),z(nz)   ! local values of x, y and z
	real(kind=8), allocatable, save        :: VALS(:,:,:,:,:)     ! (y,z,t,id,iside)
	real(kind=8), allocatable, save        :: DERIVS(:,:,:,:,:)   ! (y,z,t,id,iside)	
	real(kind=8), allocatable, save        :: TIME(:)             ! time in boundary data files
	integer                                :: iside,i0,i1
	real(kind=8)                           :: t_left,t_right,xx
	logical, save                          :: do_bottom,do_top,do_read
	
	!------------------------------------------------------------
 	! filename stuff...
 	!------------------------------------------------------------
	integer                         ::  ncid,varid,i_proc,j_proc,dimid
	character(len=4)                ::  cid
	character(len=3), save          ::  iproc,jproc
	character(len=80)               ::  ncfile,val_name,deriv_name,char_side
	integer, save                   ::  start(3),count(3),nvars=4,read_time=.TRUE.
	logical, save                   ::  first_entry=.TRUE.
	integer                         ::  i,j,icount,remaining
	include 'netcdf.inc'
	
	if( first_entry ) then
	
		if( do_second_scalar ) nvars = 5
		
		t_end = 0.d0
		start_slice = 1
		
		allocate(   VALS(nx,ny,slices_per_chunk,nvars,2) )  ! B/T (x,y,t,id,side)
		allocate( DERIVS(nx,ny,slices_per_chunk,nvars,2) )
		VALS = 0.d0 ; DERIVS = 0.d0	
		allocate( TIME(slices_per_chunk) )
		TIME = 0.d0	
		
		write(unit=iproc,fmt='(I3.3)') proc_row(YBLOCK,myid) ! 3 digit char string for iproc
   		write(unit=jproc,fmt='(I3.3)') proc_col(YBLOCK,myid) ! 3 digit char string for jproc
   		
   		!-----------------------------------------------------------
		! determine if this processor needs to fill tha data values
		!-----------------------------------------------------------
		do_bottom=.FALSE. ;  do_top=.FALSE.
		if( z(1) == 0.d0 )   do_bottom = .TRUE.
		if( z(nz) == Lz )    do_top = .TRUE.
		
		if( myid == 0 ) then
			open(1,file='output/debug_data/BC_data_BT_specs')
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) ' initial and final integration times [s]   :',t0,tf
				write(1,*) ' processor myid=0  iproc,jproc             :',iproc,' , ',jproc
				write(1,*) ' chunk size nx, ny, slices_per_chunk,nvars :',nx,ny,slices_per_chunk,nvars
				write(1,*) ' for this processor do_bottom,do_top       :',do_bottom,do_top
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) '---------------------------------------------------------------------------'
				write(1,*) ' '
			close(1)	
		endif
   	
		first_entry=.FALSE.
	endif
	   
	if( id == 1 ) then
		val_name = 'u' ; deriv_name = 'u_z'
	elseif( id == 2 ) then
		val_name = 'v' ; deriv_name = 'v_z'
	elseif( id == 3 ) then
		val_name = 'w' ; deriv_name = 'w_z'
	elseif( id == 4 ) then
		val_name = 's1' ; deriv_name = 's1_z'
	elseif( id == 5 ) then
		val_name = 's2' ; deriv_name = 's2_z'
	endif
			
	!-----------------------------------------------------------
	! determine if a new data chunk needs to be read and if so,
	! read the new chunk of boundary data
	!-----------------------------------------------------------
	if( do_bottom .or. do_top ) then
		
		if( t >= t_end ) then    ! a new chunk of boundary data is needed, read it
	
			do iside = 1,2    ! B to T
   		
   				if( iside==1 .and. do_bottom ) then
   					if(proc_col(YBLOCK,myid) .NE. 0) stop 'problem detecting do_bottom in bvals_BT' 
   					ncfile = 'BVALS/bottom_'//iproc//'-'//jproc//'.nc'  ! i.e.bottom_013-000.nc
   					do_read = .TRUE.   					
   				elseif( iside==2 .and. do_top ) then
   					if(proc_col(YBLOCK,myid) .NE. p2-1) stop 'problem detecting do_top in bvals_BT'
   					ncfile = 'BVALS/top_'//iproc//'-'//jproc//'.nc'  ! i.e. top_000-008.nc
   					do_read = .TRUE.
   				else
   					do_read = .FALSE.
   				endif
   			
   				if( do_read ) then
   					!--------------------------------------
   					!  open the file, check for error
   					!--------------------------------------
   					ierr=NF_OPEN(ncfile,NF_NOWRITE,ncid)
   					if(ierr.ne.NF_NOERR) then
    					write(0,*) '... ERROR OPENING NETCDF FILE: bvals_BT ',trim(ncfile)
    					write(0,*) '... myid, ierr ',myid, ierr
    					stop
   					endif
   			
   					!-----------------------------------------
   					!  determine total number of time slices
   					!-----------------------------------------
   					ierr = NF_INQ_DIMID(ncid, "timedimension", dimid)
					ierr = NF_INQ_DIMLEN(ncid,dimid,num_time_slices)
					if(ierr.ne.NF_NOERR) then
    					write(0,*) '... ERROR GETTING TIMEDIMENSION LENGTH: bvals_BT ',trim(ncfile)
    					write(0,*) '... myid, ierr ',myid, ierr
    					stop
   					endif
					if( num_time_slices < slices_per_chunk ) then
						slices_per_chunk = num_time_slices  ! don't read more slices than there are in the file
					endif
					
					!---------------------------------------------
   					!  set start and count to read desired data
   					!---------------------------------------------
   					remaining = num_time_slices - start_slice + 1
   					icount = min(remaining,slices_per_chunk)
					start=(/1,1,start_slice/)
   					count=(/nx,ny,icount/)
   							
					!-----------------------------------------
   					!  get start and end times for this chunk
   					!  set the time increment
   					!-----------------------------------------
					ierr = NF_INQ_VARID(ncid,'time',varid)
					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start_slice,icount,TIME)
					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR GET_VARA_DOUBLE -> time: ',ierr,start_slice,icount    					
    					stop
   					endif
					t_start = TIME(1)
					dt_bvals = TIME(2)-TIME(1)
					
					if( t < TIME(1) .or. t > TIME(icount) ) then
						if(myid==0) then
							write(0,*) id,myid," time range of new chunk does not span current time t"
							write(0,*) id,myid,"       t,TIME(1),TIME(icount): ",t,TIME(1),TIME(icount)
							write(0,*) id,myid,"       start_slice,icount,num_time_slices: ",start_slice,icount,num_time_slices
						endif	
					endif
												   			
   					!------------------------------------------------
   					!  extract variable id for VAL , check for error
   					!------------------------------------------------
   					ierr=NF_INQ_VARID(ncid,val_name,varid)
   					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR INQ_VARID -> val, id: ',id,ierr
    					stop
   					endif

   					!------------------------------------------------
   					!  read the corresponding variable (dimensional)
   					!------------------------------------------------
   					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start,count,VALS(1,1,1,id,iside)) 
   					if (ierr.ne.NF_NOERR) then
   						write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> val, id: ',id,ierr
    					stop
   					endif
   				
   					!---------------------------------------------------
   					!  extract variable id for DERIV , check for error
   					!---------------------------------------------------
   					ierr=NF_INQ_VARID(ncid,deriv_name,varid)
   					if (ierr.ne.NF_NOERR) then
    					write(0,*) myid,'NetCDF ERROR INQ_VARID -> deriv, id: ',id,ierr
    					stop
   					endif

   					!------------------------------------------------
   					!  read the corresponding variable (dimensional)
   					!------------------------------------------------
   					ierr = NF_GET_VARA_DOUBLE(ncid,varid,start,count,DERIVS(1,1,1,id,iside)) 
   					if (ierr.ne.NF_NOERR) then
   						write(0,*) myid,'NetCDF ERROR NF_GET_VARA_REAL -> deriv, id: ',id,ierr
    					stop
   					endif
   				
   					if( myid == 0 .and. id<9 ) then
   						open(1,file='output/debug_data/BC_data_BT_specs',position='append')
   						write(1,*) ' finished reading file, iside        :',trim(ncfile),iside
   						write(1,*) ' id and variable names read          :',id,trim(val_name),' , ',trim(deriv_name)
   						write(1,*) ' total number of time slices in file :',num_time_slices
   						write(1,*) ' slice range for this chunk          :',start_slice,start_slice + icount - 1 
   						write(1,*) ' time range for this chunk           :',TIME(1),TIME(icount)
   						write(1,*) ' time value triggering read          :',t
   						write(1,*) '  '
   						write(1,*) '  '
   						close(1)
   					endif
   				
   				endif ! do_read block
   			
   			enddo   ! end loop over 2 endpoints
   			   			   		   			
		endif   ! end block for reading a new chunk of boundary data
		
		!-----------------------------------------------------------------------
		! a valid boundary data chunk is now available (new or old)
		! interpolate between saved BC values to get values at the current time
		!-----------------------------------------------------------------------
		
   		i0 = floor( (t-t_start)/dt_bvals ) + 1        ! index of saved time just lt t
   		i1 = minval( (/i0 + 1,slices_per_chunk/) )    ! index of saved time just gt t
   		t_left = t_start + (i0-1)*dt_bvals            ! saved time in [s] just lt t
   		t_right = t_left + dt_bvals                   ! saved time in [s] just gt t
   		xx = (t - t_left)/dt_bvals                    ! linear interpolation weight
   		   
   		!---------------
   		! SANITY CHECK
   		!---------------
   		if( t < t_left .OR. t > t_right ) then
    		if(myid==0) then
     			write(0,*) 'misalignment in BC time values: ',t,t_start,t_end,dt_bvals
     			write(0,*) 't_left, t_right ',t_left, t_right
     			write(0,*) 'i0, i1, xx ',i0,i1,xx
     			STOP
    		endif
   		endif   
   		if( xx < 0.d0 .OR. xx > 1.d0 ) then
    		if(myid==0) then
     			write(0,*) 'misalignment in BC time values: ',t,t_start,t_end
     			write(0,*) 't, t_left, t_right ',t, t_left, t_right
     			write(0,*) 'xx, dt_bvals  ',xx, dt_bvals
     			STOP
    		endif
   		endif
   		   		   		 
   		!-------------------------------------------------
   		! do the interpolation and store the value
   		!-------------------------------------------------
   		if( do_bottom ) then
   			iside = 1
   			bottom_vals(:,:,id)   =   VALS(:,:,i0,id,iside) + xx*(   VALS(:,:,i1,id,iside) -   VALS(:,:,i0,id,iside) )  ! time interpolated values
   			bottom_derivs(:,:,id) = DERIVS(:,:,i0,id,iside) + xx*( DERIVS(:,:,i1,id,iside) - DERIVS(:,:,i0,id,iside) )  ! time interpolated derivs
   		endif
   		
   		if( do_top ) then
   			iside = 2
   			top_vals(:,:,id)   =   VALS(:,:,i0,id,iside) + xx*(   VALS(:,:,i1,id,iside) -   VALS(:,:,i0,id,iside) )  ! time interpolated values
   			top_derivs(:,:,id) = DERIVS(:,:,i0,id,iside) + xx*( DERIVS(:,:,i1,id,iside) - DERIVS(:,:,i0,id,iside) )  ! time interpolated derivs
   		endif
   		
   		
	endif   ! end block for processors that own EW boundaries
	
	!-------------------------------------------------
   	! workload can be different, not all processors
   	! have boundaries and do the reads, sync back up
   	!-------------------------------------------------		
	call mpi_barrier(comm,ierr)  
	
 return	          
end subroutine bvals_BT

