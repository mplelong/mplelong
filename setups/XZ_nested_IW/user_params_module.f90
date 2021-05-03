!!=============================================
!!   Module for user parameters and variables
!!   to access variables, routines etc
!!   ==>   use user_params
!!=============================================

module user_params

	!-------------------------------------------------------------------------------------
	! Parameters defining the "parent soln" internal wave mode
	!  (f0, g, rho0 all set in problem_params)
	!-------------------------------------------------------------------------------------
	
	real(kind=8), parameter          :: L = 1.5d5         !! [m]
	real(kind=8), parameter          :: H = 3000.d0       !! [m]
	real(kind=8), parameter          :: N = 2.d-3         !! [1/s]  used in setting ambient s1=rho_bar(z) 
	real(kind=8), parameter          :: mode_x = 1.0d0    !! periodic mode number in L
	real(kind=8), parameter          :: mode_y = 0.0d0    !! periodic mode number in L
	real(kind=8), parameter          :: mode_z = 1.0d0    !! cos mode number in H
	real(kind=8), parameter          :: x0 = 0.5d0*L      !! offset for corner of computational domain
	real(kind=8), parameter          :: y0 = 0.4d0*L      !! offset for corner of computational domain
	real(kind=8), parameter          :: z0 = 0.6d0*H      !! offset for corner of computational domain
	
	!-------------------------------------------------------------------------------------
	! Parameters defining additional interior body forcing
	!-------------------------------------------------------------------------------------
	
	real(kind=8), parameter          :: Lx=30000.d0       !! MUST MATCH VALUE IN PROBLEM_PARAMS
	real(kind=8), parameter          :: Ly=30000.d0       !! MUST MATCH VALUE IN PROBLEM_PARAMS
	real(kind=8), parameter          :: Lz=600.d0         !! MUST MATCH VALUE IN PROBLEM_PARAMS
	real(kind=8), parameter          :: f0=1.d-4          !! MUST MATCH VALUE IN PROBLEM_PARAMS


	real(kind=8), parameter          :: xf = Lx/2.d0      !! [m] center of forcing zone
	real(kind=8), parameter          :: zf = Lz/2.d0      !! [m] center of forcing zone
	real(kind=8), parameter          :: gamma_x=Lx/6.     !! [m] size of forcing zone in x
	real(kind=8), parameter          :: gamma_z=Lz/128.   !! [1] size of forcing zone in z
	real(kind=8), parameter          :: sigma = f0/5.     !! [1/s] frequency of forcing
	real(kind=8), parameter          :: A = 0.*f0         !! [m/s2] amplitude of acceleration term on rhs of u


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
		implicit none
		!--------------------------------------------------------------------------------------------
		!  For example, can set values different than defaults for Bernoulli/Cosine differentiation
		!  Default values  Q=5 (4 term expansions)  and filter_fraction = 0.05
		!  	filter_fraction = 0 ==> no spectral filtering;  filter_fraction < 0 ==> 2/3 rule
		!  This is invoked by the "call reset_default_values" in user_ICS.
		!--------------------------------------------------------------------------------------------
		Q = 9                                ! gives best results for differentiation tests                      
		filter_fraction = 0.000              ! .025 using 2/3 rule gives worse results near E/W boundaries
		
		variable_key(9:11,:) = 1             ! for debugging, output ustar,vstar and wstar
		
	end subroutine change_default_values


	real(kind=8) function parent_soln(X,Y,Z,t,id)
	!-----------------------------------------------------------------
	! function to evaluate parent soln at X,Z,t  ; 
	! here an IW mode with wave parameters available via user_params
	! id = [1,2,3,4,6] for [U,V,W,PD,P] respectively
	!-----------------------------------------------------------------
		use dimensional_scales,  only: f0,rho0,g
		implicit none
		real(kind=8), intent(in)    :: X,Y,Z,t
		integer, intent(in)         :: id
		real(kind=8),save           :: A,kx,ky,kz,omega,omega2,phase,N2,f2
		real(kind=8)                :: pi,argxy,argz,ans
		logical, save               :: first_entry=.TRUE.
		if( first_entry ) then
			pi = 4.d0*atan(1.d0)
			A = 0.05    ! 0.01 in most of the tests
			kz = mode_z * pi/H
			kx = mode_x * 2.d0*pi/L
			ky = mode_y * 2.d0*pi/L
			f2 = f0*f0 
			N2 = N*N
			omega2 = ((kx**2+ky**2)*N2 + kz**2*f2)/(kx**2 + ky**2 + kz**2)
			omega = sqrt(omega2)
			phase = pi/4.d0	
			if(ky .NE. 0.d0 )	stop 'Analytical wave mode soln only valid for ky=0 (rederive for general case)'	
			first_entry = .FALSE.
		endif
		
		argxy = kx*X + ky*Y - omega*t + phase
		argz = kz*Z
		
		if(id==1) ans = A*cos(argz)*cos(argxy)                              ! U
		if(id==2) ans = A*(f0/omega)*cos(argz)*sin(argxy)                   ! V
		if(id==3) ans = A*(kx/kz)*sin(argz)*sin(argxy)                      ! W
				
		if( id==4 ) then
			ans = -A*(kx/kz)*(N2/omega)*sin(argz)*cos(argxy)                ! B
			ans = -(rho0/g)*ans  
		endif
		
		if(id==6) ans = -A*(1./(kx*omega))*(f2-omega2)*cos(argz)*cos(argxy) ! P
	
		parent_soln = ans
	 return	
	end function parent_soln
	
	
	
	real(kind=8) function parent_deriv(X,Y,Z,t,id,dir)
	!-----------------------------------------------------------------
	! function to evaluate derivs of the parent soln at X,Z,t  ; 
	! here an IW mode with wave parameters available via user_params
	! id = [1,2,3,4,6] for [U,V,W,PD,P] respectively
	!-----------------------------------------------------------------
		use dimensional_scales,    only: f0,rho0,g
		implicit none
		real(kind=8), intent(in)      :: X,Y,Z,t
		integer, intent(in)           :: id
		character(len=80), intent(in) :: dir
		real(kind=8),save             :: A,kx,ky,kz,omega,omega2,phase,N2,f2
		real(kind=8)                  :: pi,argxy,argz,ans
		logical, save                 :: first_entry=.TRUE.
		if( first_entry ) then
			pi = 4.d0*atan(1.d0)
			A = 0.05    ! 0.01 in most of the tests
			kz = mode_z * pi/H
			kx = mode_x * 2.d0*pi/L
			ky = mode_y * 2.d0*pi/L
			f2 = f0*f0 
			N2 = N*N
			omega2 = ((kx**2+ky**2)*N2 + kz**2*f2)/(kx**2 + ky**2 + kz**2)
			omega = sqrt(omega2)
			phase = pi/4.d0
			first_entry = .FALSE.
		endif
		
		argxy = kx*X + ky*Y - omega*t + phase
		argz = kz*Z
		
		if(id==1) then         ! U = A*cos(argz)*cos(argxy)
			if(dir=='x') ans= -kx * A*cos(argz)*sin(argxy)
			if(dir=='y') ans= -ky * A*cos(argz)*sin(argxy)
			if(dir=='z') ans= -kz * A*sin(argz)*cos(argxy)
		elseif(id==2) then     ! V = A*(f0/omega)*cos(argz)*sin(argxy)
			if(dir=='x') ans=  kx * A*(f0/omega)*cos(argz)*cos(argxy)
			if(dir=='y') ans=  ky * A*(f0/omega)*cos(argz)*cos(argxy)
			if(dir=='z') ans= -kz * A*(f0/omega)*sin(argz)*sin(argxy)
		elseif(id==3) then     ! W = A*(kx/kz)*sin(argz)*sin(argxy)
			if(dir=='x') ans=  kx * A*(kx/kz)*sin(argz)*cos(argxy)
			if(dir=='y') ans=  ky * A*(kx/kz)*sin(argz)*cos(argxy)
			if(dir=='z') ans=  kz * A*(kx/kz)*cos(argz)*sin(argxy)
		elseif(id==4) then     ! B = -A*(kx/kz)*(N2/omega)*sin(argz)*cos(argxy)
			if(dir=='x') ans=  kx * A*(kx/kz)*(N2/omega)*sin(argz)*sin(argxy)
			if(dir=='y') ans=  ky * A*(kx/kz)*(N2/omega)*sin(argz)*sin(argxy)
			if(dir=='z') ans= -kz * A*(kx/kz)*(N2/omega)*cos(argz)*cos(argxy)
			ans = -(rho0/g)*ans    ! convert buoyancy B to perturbation density:  b = -(g/rho0)*pd
		endif		                               				
		parent_deriv = ans
	 return	
	end function parent_deriv
	
end module user_params

