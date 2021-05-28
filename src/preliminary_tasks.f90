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
	use methods_params
	use etc
	use user_params, only: change_default_values
	implicit none
	integer                 :: i,n_seed
	integer, allocatable    :: seed(:)
	
	
	!-------------------
	!  Initialize mpi
	!-------------------
 	call start_mpi(myid,np,logfile,comm)
 	numprocs=np
 	
 	!-----------------------------------------------------------
	! Initialize random number generator
	! F90 intrinsics RANDOM_SEED, RANDOM_NUMBER
	! --> make sure each processor gets a unique random sequence
	!-----------------------------------------------------------
	call RANDOM_SEED(size = n_seed)   ! get the size of the seed sequence
	allocate(seed(n_seed))            ! create a temp array to store the seed sequence

 	seed(1) = myid                    ! starting point unique to each processor
 	do i=2,n_seed
  		seed(i) = seed(i-1) + 1       ! some arbitrary unique sequence for each processor
 	enddo
 	call RANDOM_SEED(put=seed(1:n_seed))   ! initialize the generator with this sequence       
 	deallocate(seed)                       !  don't need this array anymore
 
 	
 	!-------------------
	!  Read user data
	!-------------------
	if(myid==0) then
  		message = '................  Reading input/user_params and input/io_params ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	call ReadUserData
 	
 	!------------------------------------------------------------
 	! Change any default values specified in user_params_module
 	!------------------------------------------------------------
  	call change_default_values

 	
 	!---------------------------------------
	!  Set params for data decomposition
	!---------------------------------------
	if(myid==0) then
  		message = '................  Calling Decomposition2D ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
 	call Decomposition2D
 	 	
 	
 	!-----------------------------------------------
	! Allocate dependent and intermediate variables
	!-----------------------------------------------	
	if(myid==0) then
  		message = '................  Allocating memory for dependent/intermediate variables ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif	
 	call AllocateIndependentVariables
 	call AllocateDependentVariables
 	call AllocateIntermediateVariables
 	
 	
 	!------------------------------------------
	! test transposes for hangs w/ zero arrays
	!------------------------------------------
 	call test_transposes
 	

 	!-------------------------------------
	! global coordinate arrays x,y and z
	!-------------------------------------
	if(myid==0) then
  		message = '................  Initializing independent variables ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
	call initialize_coord(x,nx,Lx,x_periodic)
	call initialize_coord(y,ny,Ly,y_periodic)
	call initialize_coord(z,nz,Lz,z_periodic)
		
	
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
	
	!--------------------------
	!  set some initial values
	!--------------------------
	t_secs = t0
	tn = t0
	tnp1 = tn + dt
	istep = istart
	
	!-----------------------------------------------
	! user may have called "change_default_values"
	! check for FS rigid lid xy periodic
	!-----------------------------------------------
	if( FS_XY_PERIODIC ) then
		x_periodic = .TRUE.
		y_periodic = .TRUE.
		z_periodic = .FALSE.
		z_FSRL = .TRUE.
		Q = 0
		user_bcs = .FALSE.
		endpoint_smoothing = .FALSE.
	endif
	
	!-----------------------------------------------
	! In this case, the Q=0 setting implies cos in y
	! no BCs needed at y=0,Ly (or anywhere else)
	!-----------------------------------------------
	if( z_FSRL .and. x_periodic .and. Q==0 ) then
		user_bcs = .FALSE.
		endpoint_smoothing = .FALSE.
	endif
	
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
	! get t=t0 boundary conditions and apply
	! them to the initial conditions as a test
	!-----------------------------------------------
 	tnp1 = tn	!  set to tn temporarily for testing
 	call fill_boundary_arrays
 	call apply_bcs		
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
 	call initialize_fourier_stuff
 	

	
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
	
	if(myid==0) then
  		message = '................  Verifying Bernoulli Cosine differentiation ...'
  		write(0,*) message
  		call LogMessage(message,logfile)
 	endif
	call verify_deriv_BC
		 	 	
 	!-------------------------------------------
	! test some essential higher order routines
	!-------------------------------------------
	call test_divergence
	call test_gradient
	call test_mudotgradf
 	call test_transform_xy
 	call test_poisson_solver
 
	
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
