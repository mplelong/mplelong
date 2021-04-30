!--------------------------------------------------------------------------!
!    Subroutine to do preliminary tasks to prepare for flow_solve
!    integration in time. Uses Bernoulli/Cosine method for differentiation.
!--------------------------------------------------------------------------!
subroutine preliminary_tasks
	use mpi_params
	use independent_variables
	use dependent_variables
	use intermediate_variables
	use decomposition_params
	use differentiation_params
	use fourier_differentiation_tools
	use methods_params,only             : high_order_operators
	use etc
	implicit none
	logical                            :: debug
	character(len=80)                  :: exp_type
	integer                            :: dim
	integer(kind=8)                    :: plan_i
	
	real(kind=8)                       :: pi,diff,tol
	integer(kind=8)                    :: plans(2)
	integer                            :: i,j,k,dir,order,Qval
	integer                            :: locnx, locnz
	real(kind=8),allocatable           :: in(:),out(:),myxvals(:),myzvals(:)
	
	
	pi = 4.d0*atan(1.d0)
	
	!-------------------
	!  Initialize mpi
	!-------------------
 	call start_mpi(myid,np,logfile,comm)
 	numprocs=np
 
 	
 	!-------------------
	!  Read user data
	!-------------------
	if(myid==0) then
  		message = '...  Reading input/user_params and input/io_params ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	call ReadUserData
 	
 	!--------------------------
	!  set some initial values
	!--------------------------
	tn = t0
	tnp1 = tn + dt

 	
 	!---------------------------------------
	!  Set params for data decomposition
	!---------------------------------------
	if(myid==0) then
  		message = '...  Calling Decomposition2D ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	call Decomposition2D
 	
 	locnx = array_size(JDIM,YBLOCK,myid)
	locnz = array_size(KDIM,YBLOCK,myid)
	allocate( myxvals(locnx), myzvals(locnz) )
	myxvals = 0.d0 ; myzvals = 0.d0
 	
 	
 	!-----------------------------------------------
	! Allocate dependent and intermediate variables
	!-----------------------------------------------	
	if(myid==0) then
  		message = '...  Allocating memory for dependent/intermediate variables ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	
 	call AllocateIndependentVariables
 	call AllocateDependentVariables
 	call AllocateIntermediateVariables
 	
 	
 	if(myid==0) then
  		message = ' ................  Calling transpose routines with zero arrays to test for hangs ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
  		tmpX =1.d0 ; tmpY=2.d0; tmpZ=3.d0
 	endif
 	call xblock_2_yblock(tmpX,tmpY)
 	call yblock_2_xblock(tmpY,tmpX)
 	call yblock_2_zblock(tmpY,tmpZ)   
 	call zblock_2_yblock(tmpZ,tmpY)
 	if(myid==0) then
  		message = ' ................                      X <--> Y and Y <--> Z execute  ... '
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	

 	!-------------------------------------
	! global coordinate arrays x,y and z
	!-------------------------------------
	if(myid==0) then
  		message = '................  Initializing coordinates ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
	call initialize_coord(x,nx,Lx)
	call initialize_coord(y,ny,Ly)
	call initialize_coord(z,nz,Lz)
	
	call get_my_xvals(myxvals,YBLOCK,myid)
	call get_my_zvals(myzvals,YBLOCK,myid)
	
	
	!-----------------------------------------------
	! call user routine to prescribe the initial
	! conditions for u,v,w and scalar(s)
	!-----------------------------------------------	
	if(myid==0) then
  		message = '................  Calling user routine for initial conditions ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif 	
	call InitialConditions
	
	!-----------------------------------------------
	! call user routine to prescribe the ambient
	! scalar profiles s1_bar(z) and s2_bar(z)
	!-----------------------------------------------	
	if(myid==0) then
  		message = '................  Calling user routine to prescribe scalar profiles s1_bar(z) and s2_bar(z) ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif 	
	call prescribe_ambient_scalar_profiles
	
	
	!-----------------------------------------------
	! get t=0 boundary conditions and apply
	! them to the initial conditions as a test
	!-----------------------------------------------
 	tnp1 = tn	!  set to 0 temporarily for testing
 	call fill_boundary_arrays
 	call apply_bcs
 	if(myid==0) then
 		open(1,file='output/debug_data/check_NS_bdry_scheme')
 			do j=1,ny
 				write(1,*) y(j),u(j,1,1),v(j,1,1),w(j,1,1),s1(j,1,1)
 			enddo	
 		close(1)
 	endif
 	if( myzvals(1)==0.d0 ) then
 		open(1,file='output/debug_data/check_B_bdry_scheme')			
 			do k=1,locnz
 				write(1,*) myzvals(k),u(1,1,k),v(1,1,k),w(1,1,k),s1(1,1,k)
 			enddo
 		close(1)
 	endif
 	if( myzvals(locnz)==Lz ) then	 		
 		open(1,file='output/debug_data/check_T_bdry_scheme')
 			do k=1,locnz
 				write(1,*) myzvals(k),u(1,1,k),v(1,1,k),w(1,1,k),s1(1,1,k)
 			enddo
 		close(1)
 	endif 		
 	tnp1 = tn + dt    !  reset to proper value
	
	
	!----------------------------------------------------
	! setup for cos expansions in x,y & z
	! fill wavenumber and filter arrays
	! create FFTW3 plans for cos and sin
	!---------------------------------------------------- 	
 	if(myid==0) then
  		message = '................  Initializing wavenumber/filter arrays, creating FFTW3 plans ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	
 	allocate( kx(nx),ky(ny),kz(nz) )
 	allocate( kxfilter(nx), kyfilter(ny), kzfilter(nz) )
 	
 	dim = 1   ! x
 	exp_type = 'cos'
	call fourier_init(exp_type,nx,Lx,kx,kxfilter,cos_plan(dim),plan_i)
	exp_type = 'sin'
	call fourier_init(exp_type,nx,Lx,kx,kxfilter,sin_plan(dim),plan_i)
	if(myid==0) then
		open(1,file='output/debug_data/x_wavenumbers')
		do i=1,nx
			write(1,*) i,kx(i),kxfilter(i)
		enddo
		close(1)
	endif
	
	dim = 2   ! y
 	exp_type = 'cos'
	call fourier_init(exp_type,ny,Ly,ky,kyfilter,cos_plan(dim),plan_i)
	exp_type = 'sin'
	call fourier_init(exp_type,ny,Ly,ky,kyfilter,sin_plan(dim),plan_i)
	if(myid==0) then
		open(1,file='output/debug_data/y_wavenumbers')
		do j=1,ny
			write(1,*) j,ky(j),kyfilter(j)
		enddo
		close(1)
	endif
	
	dim = 3   ! z
 	exp_type = 'cos'
	call fourier_init(exp_type,nz,Lz,kz,kzfilter,cos_plan(dim),plan_i)
	exp_type = 'sin'
	call fourier_init(exp_type,nz,Lz,kz,kzfilter,sin_plan(dim),plan_i)
	if(myid==0) then
		open(1,file='output/debug_data/z_wavenumbers')
		do k=1,nz
			write(1,*) k,kz(k),kzfilter(k)
		enddo
		close(1)
	endif
	
	!----------------------------------------------------
	! verify fourier cosine differentiation
	!----------------------------------------------------
	tol = 1.d-10
	allocate( in(nx),out(nx) )
		in=0.d0 ; out=0.d0
 		do j=1,nx
 			in(j) = cos( 2.*pi*x(j)/Lx )
 		enddo
 	
 		dir = 1
 		plans(1) = cos_plan(dir)
 		plans(2) = sin_plan(dir)
 		exp_type = 'cos'
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
 		plans(1) = cos_plan(dir)
 		plans(2) = sin_plan(dir)
 		exp_type = 'cos'
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
 		plans(1) = cos_plan(dir)
 		plans(2) = sin_plan(dir)
 		exp_type = 'cos'
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
 	 	
	
 	if(myid==0) then
  		message = '................  Calculating diffusion coefficients ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	call diffusion_coeffs
	
	
	!----------------------------------------------------
	! evaluate basis functions for
	! the Bernoulli expansions in x,y & z
	! results available in module differentiation_params
	!----------------------------------------------------
	if(myid==0) then
  		message = '................  Evaluating/storing Bernoulli basis functions ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
	call evaluate_basis_functions
	
	
	
		
	!----------------------------------------------------
	! fill and factor the coefficient matices for
	! computing the B. expansion coefficients in x,y & z
	! results available in module differentiation_params
	!----------------------------------------------------
	if(myid==0) then
  		message = '................  Forming/factoring matrices for Bernoulli expansions ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
	call setup_factor_B_matrices
	

	
	!----------------------------------------------------
	! verify BC differentiation
	!----------------------------------------------------
	allocate( in(nx),out(nx) )
		debug = .FALSE.
		in=0.d0 ; out=0.d0
 		do j=1,nx
 			in(j) = sin( pi*x(j)/Lx )   ! even extension has discontinuous derivs at x=0, x=Lx
 		enddo
 	
 		dir = 1
 		Qval = Q
 		call deriv_BC(in,out,nx,dir,Qval,debug)  ! Bernoulli/Cosine derivative
 	
 		if(myid==0 .and. nx > 1) then
 			do j=1,nx
 				diff = abs(out(j) - (pi/Lx)*cos( pi*x(j)/Lx ))
 				!write(0,*) x(j), in(j), out(j), (pi/Lx)*cos( pi*x(j)/Lx )
 				if( diff  > tol ) stop ' problem verifying x Bernoulli-cosine differentiation in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) '                         ....................... d/dx using Bernoulli-cosine method looks fine, tol=',tol
 	deallocate( in,out )
 	
	allocate( in(ny),out(ny) )
		debug = .FALSE.
		in=0.d0 ; out=0.d0
 		do j=1,ny
 			in(j) = sin( pi*y(j)/Ly )   ! even extension has discontinuous derivs at y=0, y=Ly
 		enddo
 	
 		dir = 2
 		Qval = Q
 		call deriv_BC(in,out,ny,dir,Qval,debug)  ! Bernoulli/Cosine derivative
 	
 		if(myid==0 .and. ny> 1) then
 			do j=1,ny
 				diff = abs(out(j) - (pi/Ly)*cos( pi*y(j)/Ly ))
 				!write(0,*) y(j), in(j), out(j), (pi/Ly)*cos( pi*y(j)/Ly )
 				if( abs(diff) > tol ) stop ' problem verifying y Bernoulli-cosine differentiation in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) '                         ....................... d/dy using Bernoulli-cosine method looks fine, tol=',tol
 	deallocate( in,out )
 	
 	allocate( in(nz),out(nz) )
 		debug = .TRUE.
		in=0.d0 ; out=0.d0
 		do j=1,nz
 			in(j) = sin( pi*z(j)/Lz )   ! even extension has discontinuous derivs at z=0, z=Lz
 		enddo
 	
 		dir = 3
 		Qval = Q
 		call deriv_BC(in,out,nz,dir,Qval,debug)  ! Bernoulli/Cosine derivative
 	
 		if(myid==0 .and. nz>1) then
 			do j=1,nz
 				diff = abs(out(j) - (pi/Lz)*cos( pi*z(j)/Lz ))
 				!write(0,*) z(j), in(j), out(j), (pi/Lz)*cos( pi*z(j)/Lz ), diff
 				if( diff > tol ) stop ' problem verifying z Bernoulli-cosine differentiation in preliminary_tasks '
 			enddo
 		endif
 		write(0,*) '                         ....................... d/dz using Bernoulli-cosine method looks fine, tol=',tol
 	deallocate( in,out ) 	 	
 	
 	
 	!-------------------------------------------
	! test some essential higher order routines
	!-------------------------------------------
	call test_divergence
	call test_gradient  
	call test_udotgradf
 	call test_transform_xy
 	call test_poisson_solver
 	call test_z_diffusion
 
	
	if(myid==0) then
		message = '................   -----> routine PreliminaryTasks exiting normally <----------  ................ '
		write(0,*) ' '
		write(0,*) '-------------------------------------------------------------------------------------------------'
		write(0,*) message
		write(0,*) '-------------------------------------------------------------------------------------------------'
		write(0,*) ' '
		call LogMessage(message,logfile)
	endif 

	call mpi_barrier(comm,ierr)
 return	
end subroutine preliminary_tasks
