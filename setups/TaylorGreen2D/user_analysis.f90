subroutine user_analysis
	use mpi_params,             only: myid
	use decomposition_params
	use dependent_variables,    only: u,v,w
	use methods_params,         only: do_second_scalar
	use independent_variables,  only: x,y,z,t_secs,tf,t0,dt
	use etc,                    only: istep,istart,iend
	use user_params
  
	implicit none
	integer,save                    ::  i,j,k,n,nsteps,ig,kg,id
	integer,save                    ::  locnx,locny,locnz
	real(kind=8),allocatable,save   ::  uvw(:,:),exact(:,:),t(:)
	real(kind=8), save              ::  xval,yval,zval
	logical                         ::  first_entry=.TRUE. 
	include 'mpif.h'
 
 	
	if( array_size(JDIM,YBLOCK,myid) == 1 ) return  ! no need to x avg if nx=1
  
	if( first_entry ) then
		locnx = array_size(JDIM,YBLOCK,myid)
		locny = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		
		! "randomly chosen" output location
		i = int(locnx/3)
		j = int(locny/6)
		k = int(locnz/2)
		
		ig = global_x_indices(START,YBLOCK,myid) + i - 1  ! global i index
		kg = global_z_indices(START,YBLOCK,myid) + k - 1  ! global k index
		
		xval = x(ig) + x0
		yval = y(j)  + y0
		zval = z(kg) + z0
		
		iend = (tf-t0)/dt + 1
		nsteps = iend - istart + 1
		allocate( uvw(nsteps,3), exact(nsteps,3), t(nsteps) )
		uvw = 0.d0
		exact = 0.d0
		t = 0.d0
		
		first_entry = .FALSE.
	endif
	
	!------------------------------
	!  add a new result to arrays
	!------------------------------
	n = istep - istart + 1
	uvw(n,1) = u(j,i,k)
	uvw(n,2) = v(j,i,k)
	uvw(n,3) = w(j,i,k)
	t(n) = t_secs
	
	!------------------------------
	! exact solutions
	!------------------------------
	do id=1,3
		exact(n,id) = parent_soln(xval,yval,zval,t(n),id)
	enddo
	
	!-----------------------------------
	! if last time step, write results
	!-----------------------------------
	if( myid==0 .and. istep==iend ) then
		open(1,file='output/TG_results')
			do n=1,nsteps
				write(1,101) xval,yval,zval,t(n),uvw(n,1),uvw(n,2),uvw(n,3),exact(n,1),exact(n,2),exact(n,3)
			enddo
		close(1)
	endif
101 format(10f14.8)  
  
 return  
end subroutine user_analysis
