!!=====================================================
!!   Module for user parameters and variables
!!   to access variables, routines etc
!!   ==>   use user_params
!!
!!  Tri-periodic vortex soln Antuono, M. (2020), JFM
!!=====================================================

module user_params

	!-------------------------------------------------------------------------------------
	! Parameters defining the "parent soln" internal wave mode
	!  (f0, g, rho0 all set in problem_params)
	!-------------------------------------------------------------------------------------
	
	real(kind=8), parameter          :: L  = 1.d0           !! [m]
	real(kind=8), parameter          :: U0 = 1.d0           !! [m/s]
	real(kind=8), parameter          :: psi = 0.d0          !! [1] translation parameter
	
		
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
		
		variable_key(6,:)=1                  ! for debugging, output div_u*
		variable_key(9:11,:) = 1             ! for debugging, output ustar,vstar and wstar
				
	end subroutine change_default_values


	real(kind=8) function exact_soln(x,y,z,t,id)
	!-----------------------------------------------------------------
	! function to evaluate parent soln at x,y,z,t  ; 
	! generally, id = [1,2,3,4,5,6] for [U,V,W,PD,S2,P] respectively
	!-----------------------------------------------------------------
		use dimensional_scales,  only: nu
		implicit none
		real(kind=8), intent(in)    :: x,y,z,t
		real(kind=8)                :: kxp,kyp,kzp    ! k * primed values
		real(kind=8), save          :: s1,s2          ! phase shifts
		integer, intent(in)         :: id
		real(kind=8),save           :: k
		real(kind=8)                :: pi,C,F,ans
		real(kind=8), external      :: myexp
		logical, save               :: first_entry=.TRUE.
		if( first_entry ) then
			pi = 4.d0*atan(1.d0)
			k = 2.d0*pi/L
			C = (4.d0/3.d0)*sqrt(2.d0/3.d0)
			s1 = (5.d0/6.d0)*pi
			s2 = (1.d0/6.d0)*pi
			first_entry = .FALSE.
		endif
		
		F = U0*C*myexp(-3.d0*nu*k**2*t)
		kxp = k*x + psi
		kyp = k*y + psi
		kzp = k*z + psi
		
		select case(id)
	
    	case(1)        !! u
    		ans = ( sin(kxp-s1)*cos(kyp-s2)*sin(kzp)           \
    		       -cos(kzp-s1)*sin(kxp-s2)*sin(kyp) )*F    
		case(2)        !! v
			ans = ( sin(kyp-s1)*cos(kzp-s2)*sin(kxp)           \
    		       -cos(kxp-s1)*sin(kyp-s2)*sin(kzp) )*F 
		case(3)        !! w
			ans = ( sin(kzp-s1)*cos(kxp-s2)*sin(kyp)           \
    		       -cos(kyp-s1)*sin(kzp-s2)*sin(kxp) )*F
    	case default
      		ans = 0.d0
		end select
		
		exact_soln = ans
	 return	
	end function exact_soln
	
		
end module user_params

