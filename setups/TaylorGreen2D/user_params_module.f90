!!=============================================
!!   Module for user parameters and variables
!!   to access variables, routines etc
!!   ==>   use user_params
!!
!!    TaylorGreen2D problem
!!=============================================

module user_params

	!-------------------------------------------------------------------------------------
	! Parameters defining the "parent soln" internal wave mode
	!  (f0, g, rho0 all set in problem_params)
	!-------------------------------------------------------------------------------------
	
	real(kind=8), parameter          :: L = 1.d0            !! [m]
	real(kind=8), parameter          :: vortex_mode = 1     !! [1] dimensionless, isotropic mode
	character(len=80), parameter     :: orientation = 'XY'  !! 'XY', 'XZ' or 'YZ' 
	
	
	!-------------------------------------------------------------------------------------
	! Parameters defining embedded subdomain
	!-------------------------------------------------------------------------------------
	real(kind=8), parameter          :: x0 = 0.0d0*L      	!! offset for corner of computational domain
	real(kind=8), parameter          :: y0 = 0.0d0*L      	!! offset for corner of computational domain
	real(kind=8), parameter          :: z0 = 0.0d0*L      	!! offset for corner of computational domain
	

	!--------------------------------------------------------------------------------------
	! restart parameters	 (only used if restart flag is .TRUE. in problem_params) 
	!                         e.g. 'RESTART/restart_'    'RESTART/XYZ_087900_'
	!--------------------------------------------------------------------------------------
	character(len=80),parameter      :: rs_basename = 'RESTART/restart_'      !! relative to $BASE directory
	logical,parameter                :: subtract_s1_bar=.FALSE.               !! TRUE if s1_bar is added to s1 in restart files
	logical,parameter                :: subtract_s2_bar=.FALSE.               !! TRUE if s2_bar is added to s2 in restart files
	logical,parameter                :: reset_dye=.FALSE.                     !! TRUE if s2 distribution reset in restart files


	!-------------------------------------------------------------------------------------
	! declare an allocatable array for storing x-averaged quantities on YZ planes
	!-------------------------------------------------------------------------------------
	real(kind=8), allocatable        :: X_means(:,:,:)
	real(kind=8), allocatable        :: rho_bar(:)           ! for storing xy average

contains

	subroutine change_default_values
		use differentiation_params,   only: Q, filter_fraction
		use io_params,                only: variable_key
		use independent_variables,    only: FS_XY_PERIODIC
		implicit none
		!--------------------------------------------------------------------------------------------
		!  For example, can set values different than defaults for Bernoulli/Cosine differentiation
		!  Default values  Q=5 (4 term expansions)  and filter_fraction = 0.05
		!  	filter_fraction = 0 ==> no spectral filtering;  filter_fraction < 0 ==> 2/3 rule
		!  This is invoked by the "call reset_default_values" in user_ICS.
		!--------------------------------------------------------------------------------------------
		Q = 9                                ! gives best results for differentiation tests                      
		filter_fraction = 0.000              ! .025 using 2/3 rule gives worse results near E/W boundaries
		
		variable_key(6,:)=1                  ! for debugging, output div_u*
		variable_key(9:11,:) = 1             ! for debugging, output ustar,vstar and wstar
		
		FS_XY_PERIODIC=.TRUE.                ! FS rigid lids, XY periodic BCs
		
	end subroutine change_default_values


	real(kind=8) function parent_soln(X,Y,Z,t,id)
	!-----------------------------------------------------------------
	! function to evaluate parent soln at X,Z,t  ; 
	! here an IW mode with wave parameters available via user_params
	! id = [1,2,3,4,5,6] for [U,V,W,PD,S2,P] respectively
	!-----------------------------------------------------------------
		use dimensional_scales,  only: nu
		implicit none
		real(kind=8), intent(in)    :: X,Y,Z,t
		integer, intent(in)         :: id
		real(kind=8),save           :: k
		real(kind=8)                :: pi,F,ans
		real(kind=8), external      :: myexp
		logical, save               :: first_entry=.TRUE.
		if( first_entry ) then
			pi = 4.d0*atan(1.d0)
			k = vortex_mode * 2.d0*pi/L	
			first_entry = .FALSE.
		endif
		
		F = myexp(-2.d0*nu*k**2*t)
		
		if( orientation == 'XY' ) then
			if(id==1) ans =  cos(k*X)*sin(k*Y)*F    ! u
			if(id==2) ans = -sin(k*X)*cos(k*Y)*F    ! v
			if(id==3) ans =  0.d0                   ! w
			if(id==4) ans =  cos(k*X)*sin(k*Y)*F    ! passive s1
		elseif( orientation == 'XZ' ) then
			if(id==1) ans =  cos(k*X)*sin(k*Z)*F    ! u
			if(id==2) ans =  0.d0                   ! v
			if(id==3) ans = -sin(k*X)*cos(k*Z)*F    ! w
			if(id==4) ans =  cos(k*X)*sin(k*Z)*F    ! passive s1
		elseif( orientation == 'YZ' ) then
			if(id==1) ans =  0.d0                   ! u
			if(id==2) ans =  cos(k*Y)*sin(k*Z)*F    ! v
			if(id==3) ans = -sin(k*Y)*cos(k*Z)*F    ! w
			if(id==4) ans =  cos(k*Y)*sin(k*Z)*F    ! passive s1
		endif
		
		parent_soln = ans
	 return	
	end function parent_soln
	
	
	
	real(kind=8) function parent_deriv(X,Y,Z,t,id,dir)
	!-----------------------------------------------------------------
	! function to evaluate derivs of the parent soln at X,Z,t  ; 
	! here an IW mode with wave parameters available via user_params
	! id = [1,2,3,4,6] for [U,V,W,PD,P] respectively
	!-----------------------------------------------------------------
		use dimensional_scales,    only: nu
		implicit none
		real(kind=8), intent(in)      :: X,Y,Z,t
		integer, intent(in)           :: id
		character(len=80),intent(in)  :: dir
		real(kind=8),save             :: k
		real(kind=8)                  :: pi,F,ans
		real(kind=8), external        :: myexp
		logical, save                 :: first_entry=.TRUE.
		if( first_entry ) then
			pi = 4.d0*atan(1.d0)
			k = vortex_mode * 2.d0*pi/L	
			first_entry = .FALSE.
		endif
		
		F = myexp(-2.d0*k**2*nu*t)
		
		if( orientation == 'XY' ) then
			if(id==1) then         ! u =  cos(k*X)*sin(k*Y)*F
				if(dir=='x') ans = -k*sin(k*X)*sin(k*Y)*F
				if(dir=='y') ans =  k*cos(k*X)*cos(k*Y)*F
				if(dir=='z') ans =  0.d0
			elseif(id==2) then     ! v = -sin(k*X)*cos(k*Y)*F
				if(dir=='x') ans = -k*cos(k*X)*cos(k*Y)*F
				if(dir=='y') ans =  k*sin(k*X)*sin(k*Y)*F
				if(dir=='z') ans =  0.d0
			elseif(id==3) then     ! w = 0
				ans = 0.d0
			elseif(id==4) then     ! s1 = cos(k*X)*sin(k*Y)*F
				if(dir=='x') ans = -k*sin(k*X)*sin(k*Y)*F
				if(dir=='y') ans =  k*cos(k*X)*cos(k*Y)*F
				if(dir=='z') ans =  0.d0
			endif
		elseif( orientation == 'XZ' ) then
			if(id==1) then         ! u =  cos(k*X)*sin(k*Z)*F
				if(dir=='x') ans = -k*sin(k*X)*sin(k*Z)*F
				if(dir=='z') ans =  k*cos(k*X)*cos(k*Z)*F
				if(dir=='y') ans =  0.d0
			elseif(id==2) then     ! v = 0
				ans = 0.d0
			elseif(id==3) then     ! w = -sin(k*X)*cos(k*Z)*F
				if(dir=='x') ans = -k*cos(k*X)*cos(k*Z)*F
				if(dir=='z') ans =  k*sin(k*X)*sin(k*Z)*F
				if(dir=='y') ans =  0.d0
			elseif(id==4) then     ! s1 = cos(k*X)*sin(k*Z)*F
				if(dir=='x') ans = -k*sin(k*X)*sin(k*Z)*F
				if(dir=='z') ans =  k*cos(k*X)*cos(k*Z)*F
				if(dir=='y') ans =  0.d0
			endif
		elseif( orientation == 'YZ' ) then
			if(id==1) then     ! u = 0
				ans = 0.d0
			elseif(id==2) then         ! v =  cos(k*Y)*sin(k*Z)*F
				if(dir=='y') ans = -k*sin(k*Y)*sin(k*Z)*F
				if(dir=='z') ans =  k*cos(k*Y)*cos(k*Z)*F
				if(dir=='x') ans =  0.d0
			elseif(id==3) then     ! w = -sin(k*Y)*cos(k*Z)*F
				if(dir=='y') ans = -k*cos(k*Y)*cos(k*Z)*F
				if(dir=='z') ans =  k*sin(k*Y)*sin(k*Z)*F
				if(dir=='x') ans =  0.d0
			elseif(id==4) then     ! s1 = cos(k*Y)*sin(k*Z)*F
				if(dir=='y') ans = -k*sin(k*Y)*sin(k*Z)*F
				if(dir=='z') ans =  k*cos(k*Y)*cos(k*Z)*F
				if(dir=='x') ans =  0.d0
			endif
		endif	
		
		parent_deriv = ans
	 return	
	end function parent_deriv
	
end module user_params

