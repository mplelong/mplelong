!!=============================================
!!   Module for user parameters and variables
!!   to access variables, routines etc
!!   ==>   use user_params
!!=============================================

module user_params

	!------------------------------------------------------------------------------------------------------------------------
	!  Parameters used to evaluate body forcing terms
	!------------------------------------------------------------------------------------------------------------------------
	real(kind=8),parameter            :: pi = 3.141592653589793d0
	real(kind=8),parameter            :: Omega = 2.d0*pi/(24.*3600.d0)    !! earth rotation frequency
	real(kind=8),parameter            :: R_earth = 6371.d3                !! earth radius [m]  ~6371 km
	real(kind=8),parameter            :: lat = 55.d0*2.d0*pi/360.d0       !! 55 N  [radians]
	logical,parameter                 :: beta_plane = .TRUE.              !! whether to add beta y to f0
	logical,parameter                 :: NT_terms = .TRUE.                !! whether to include NT dynamics or not (false to true, will adjust)


	!------------------------------------------------------------------------------------------------------------------------
	!  Parameters used to define the north and south density profiles
	!------------------------------------------------------------------------------------------------------------------------
	real(kind=8),parameter           :: sponge_scale = 20.d0             !! [km]      scale of N/S relaxation zones  
	real(kind=8),parameter           :: rho_deep = 1027.d0               !! [kg/m3]   common abyssal density
	real(kind=8),parameter           :: delta_rho_north = 2.1d0          !! [kg/m3]   bottom to top density difference, north profile
	real(kind=8),parameter           :: delta_rho_surf = 2.2d0           !! [kg/m3]   surface north to south density difference   
	real(kind=8),parameter           :: lambda_north = 1300.d0           !! [m]       exp scale height, north profile
	real(kind=8),parameter           :: lambda_star = 350.d0             !! [m]       exp depth scale for NS density difference
	real(kind=8),parameter           :: gamma_ns = 1.d0/16.d0            !! [1]       dless by Ly, ambient transition scale
	real(kind=8),parameter           :: xnoise = 1.d-3                   !! [1]       noise factor, dless by delta_rho_north
	real(kind=8),allocatable         :: rho_north(:), rho_south(:)       !! [kg/m3]   arrays to store N/S density profiles
	real(kind=8),allocatable         :: y_relax_window(:,:)              !! [1]       window defining relaxation zone near N/S ends
	real(kind=8),parameter           :: t_damp = 3600.d0                 !! [s] N/S   relaxation time scale   1800
	
	!------------------------------------------------------------------------------------------------------------------------
	!  Parameters used to define the surface wind stress and it's penetration
	!  Laplacian viscosity and diffusivity decay away from surface with scale (Lz-z_MLB)
	!------------------------------------------------------------------------------------------------------------------------
	logical,parameter                :: wind_stress = .TRUE.             !! whether to include wind stress in forcing terms
	logical,parameter                :: read_wind_data = .FALSE.         !! whether to read existing data or generate new random data
	real(kind=8),allocatable         :: z_MLB(:,:)                       !! [m]  height of mixed layer base as fn of (x,y),  YBLOCK
	real(kind=8),parameter           :: depth_MLB_initial=50.d0          !! [m]  initial mixed layer depth
	logical,parameter                :: MLB_constant=.TRUE.              !! if true, z_MLB is constant in time
	integer,parameter                :: dt_z_MLB=3600.d0                 !! [s]  time increment for recomputing z_MLB if necessary 
   
	real(kind=8),parameter           :: t_offset = -9999999.999          !! offset time for saved wind data to line up with start time t0
	real(kind=8),parameter           :: dt_rand = 60.d0                  !! [s] sampling interval used to compute random wind sequences
	real(kind=8),allocatable         :: u10(:)                           !! [m/s] random wind speeds at 10 m
	real(kind=8),allocatable         :: theta(:)                         !! [m/s] random wind directions at 10 m 
	real(kind=8),allocatable         :: tau_window(:)                    !! [1] window defining near surface forcing region
	real(kind=8),parameter           :: u10_mean = 12.d0                 !! [m/s] nominal mean wind speed at 10 m height
	real(kind=8),parameter           :: u10_stdev = 2.d0                 !! [m/s] nominal st dev of wind speed at 10 m height
	real(kind=8),parameter           :: u10_tau = 3.*3600.               !! [s] time scale of relaxation toward mean in Ornstein-Uhlenbeck process      
	real(kind=8),parameter           :: v_to_h = 1.d-4                   !! [1] factor to multiply horiz eddy vis/dif in ML to get vertical value

	!------------------------------------------------------------------------------------------------------------------------
	!  Parameters used to define the simple quadratic bottom drag
	!------------------------------------------------------------------------------------------------------------------------
	logical,parameter                :: bottom_drag = .TRUE.             !! whether to include bottom drag in forcing terms
	real(kind=8),allocatable         :: drag_window(:)                    !! [1]  window confining drag law to near-bottom zone
	real(kind=8),parameter           :: delta_bottom_drag = 50.d0        !! [m]  decay scale for Gaussian window
	real(kind=8),parameter           :: CD = 1.d-3                       !! [1]  d'less drag coefficient


	!-------------------------------------------------------------------------------------
	! declare an allocatable array for storing x-averaged quantities on YZ planes
	!-------------------------------------------------------------------------------------
	real(kind=8), allocatable        :: X_means(:,:,:)
	real(kind=8), allocatable        :: rho_bar(:)           ! for storing xy average

contains

	subroutine change_default_values
		use mpi_params,               only: myid
		use differentiation_params,   only: Q, filter_fraction
		use io_params,                only: variable_key
		use methods_params,           only: rs_basename,subtract_s1_bar,restart,add_restart_time
		use etc,                      only: istart
		use independent_variables,    only: z_FSRL,s1_z_BC
		use decomposition_params
		implicit none
		
		!----------------------------------------
		!  full depth, free-slip rigid lid ocean
		!----------------------------------------
		z_FSRL = .TRUE.
		
		!----------------------------------------
		! use cosine expansion for s1 = rho
		!----------------------------------------
		s1_z_BC = 'HOMOGENEOUS_NEUMANN'
		
		
		!-------------------------------------------------------------
		! treat y as purely cosine (x via periodicity, z via z_FSRL)
		!-------------------------------------------------------------
		filter_fraction = 0.d0 
		Q = 0                    
		
				
		!--------------------------------------------------------------------------------------
		! restart parameters	 (only set if restart flag is .TRUE. in problem_params) 
		!--------------------------------------------------------------------------------------
		if( restart ) then
			istart = 0                             ! to start time step counter at some value > 0 for continuity of datafile names
			!rs_basename = 'RESTART/XYZ_000768_'   ! relative to $FLOW_SOLVE_ROOT directory, note trailing underscore
			rs_basename  = 'RESTART/restart_'      ! default name when interpolated to new resolution
			subtract_s1_bar=.FALSE.                ! TRUE if s1_bar is added to s1 in restart files
			add_restart_time=.FALSE.               ! add time stored in restart file to t0,tf specified in problem_params
		endif
		
		variable_key(6,:) = 1                ! for debugging, output div_u*
		variable_key(9:11,:) = 1             ! for debugging, output ustar,vstar and wstar
							
	end subroutine change_default_values
	
	subroutine NS_profiles(z,rho_north,rho_south)
 		!--------------------------------------------------------
 		!  Routine to prescribe the N and S density profiles
 		!--------------------------------------------------------
 		use independent_variables,     only: Lz
 		implicit none
 		real(kind=8)                      :: z,rho_north,rho_south,d
  
 		d = Lz-z     ! depth [m]
 
 		!------------------------------------------------------------
 		!  prescribe the north profile
 		!------------------------------------------------------------
 		rho_north = rho_deep - delta_rho_north*exp(-d/lambda_north)
 
 		!------------------------------------------------------------
 		!  prescribe the south profile
 		!  ==> add enhanced near surface buoyancy
 		!------------------------------------------------------------
 		rho_south = rho_north - delta_rho_surf*exp(-d/lambda_star) 
  
	end subroutine NS_profiles


	
end module user_params

