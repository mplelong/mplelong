!!=============================================
!!   Module for user parameters and variables
!!   to access variables, routines etc
!!   ==>   use user_params
!!
!!    UBC_Stability asymmetric shear problem
!!=============================================

module user_params

	!------------------------------------------------------------------------------------------------------------------------
	!  Parameters used to define background profile
	!------------------------------------------------------------------------------------------------------------------------
   	real(kind=8),parameter        :: delta_U = 0.05d0      !! velocity difference [m/s] 
   	real(kind=8),parameter        :: delta_rho = 0.50d0    !! density difference [kg/m3] 
   	real(kind=8),parameter        :: h_rho = 0.02d0        !! thickness of density interface [m]
  	real(kind=8),parameter        :: h_u = 0.10d0          !! thickness of velocity interface [m] 
  	real(kind=8)                  :: z0_u                  !! height of max shear  [m]  (will be set to Lz/2)
  	real(kind=8)                  :: z0_rho                !! height of max density gradient  [m]  (will be set z0_u-z_off)
  	real(kind=8),parameter        :: z_off = 0.025d0       !! vertical offset of density interface from z0_u  [m]

	!------------------------------------------------------------------------------------------------------------------------
	!  Amplitude parameter for adding 3d "noise" to ics
	!------------------------------------------------------------------------------------------------------------------------
 	  real(kind=8),parameter       :: amp_pert = 1.d-2     !! ampl of density perturbations in ics [kg/m3]

	!------------------------------------------------------------------------------------------------------------------------
	! declare an allocatable array for holding rho_bar(z)     
	!------------------------------------------------------------------------------------------------------------------------
	real(kind=8), allocatable         :: rho_bar(:)        !! to hold local values of ambient density profile [kg/m3]


contains

	subroutine change_default_values
		use differentiation_params,   only: Q, filter_fraction
		use io_params,                only: variable_key
		use independent_variables,    only: FS_XY_PERIODIC,s1_z_BC,z_FSRL
		implicit none
		!--------------------------------------------------------------------------------------------
		!  For example, can set values different than defaults for Bernoulli/Cosine differentiation
		!  Default values  Q=5 (4 term expansions)  and filter_fraction = 0.05
		!  	filter_fraction = 0 ==> no spectral filtering;  filter_fraction < 0 ==> 2/3 rule
		!  This is invoked by the "call reset_default_values" in user_ICS.
		!--------------------------------------------------------------------------------------------
		filter_fraction = 0.05               ! .05 is default value
		
		variable_key(6,:)=1                  ! for debugging, output div_u*
		variable_key(9:11,:) = 1             ! for debugging, output ustar,vstar and wstar
		
		FS_XY_PERIODIC = .TRUE.              ! FS rigid lids, XY periodic BCs
		z_FSRL = .TRUE.                      ! set this here as well
		s1_z_BC = 'HOMOGENEOUS_NEUMANN'      ! default is Dirichlet
	 return	
	end subroutine change_default_values
	
	real(kind=8) function eval_rho_bar(z) 
 		!--------------------------------------------------------------------------------
 		!  routine will be called in user_ics
 		!      z          vertical position in [m]
 		!      rho_bar    [kg/m3]
 		!--------------------------------------------------------------------------------
 		use independent_variables, only: Lz
 		implicit none 
 		real(kind=8)             :: z,ans,zz,z0
 		real(kind=8), external   :: mytanh
 
 		zz = Lz/2.d0                !! shear interface location [m]
 		z0 = zz - z_off             !! density interface position, displaced downward by z_off                                                            
 		ans = (delta_rho/2.d0)*(1.d0 - mytanh(2.d0*(z - z0_rho)/h_rho))   
   
		eval_rho_bar = ans   
 	 return
	end function eval_rho_bar
	
end module user_params

