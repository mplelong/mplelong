!-----------------------------------------------------------
!   user configurable routines that need to be compiled
!-----------------------------------------------------------

subroutine prescribe_ambient_scalar_profiles
!----------------------------------------------------------------
!  Code is written to solve the eqn of motion for s1'
!  where s1 = s1_bar(z) + s1'(x,y,z,t). To do this efficiently,
!  we need d/dz and d2/dz2 of s1_bar as functions of z.
!  Setting all three to zero is perfectly reasonable.
!
!  User is responsible for getting the dimensions correct.
!  e.g. for s1 used as density, s1_bar should be in [kg/m3] etc
!       for s1 used as a passive scalar, s1_bar is dimensionless concentration
!
! output
! s1_bar, s1_bar_z, s1_bar_zz  [dimensional units]
! scalar concentration and z derivs  
!
! Similar format for the second scalar s2 if there is one.
! 
! Be sure to set the ambient_profile flags, code is slightly more efficient
! if they are known to be set to zero and can be ignored.
!----------------------------------------------------------------
 	use independent_variables,     only: z,nz,Lz              ! z is global nz array in [m]
 	use dependent_variables,       only: s1_bar,s2_bar        ! global, i.e. ===>(nz,3)
 	use methods_params,            only: do_second_scalar, ambient_profile
 	use dimensional_scales,        only: rho0,g
 	implicit none
 	integer                           :: k           ! loop index, can go from 1,nz
 	real(kind=8), parameter           :: N = 0.d0
 

 	!-------------------------------------------------------------------------
 	! For this problem, expand around linear density profile with specified N
 	!  N2 = -(g/rho0)* d/dz rho_bar  ==> d/dz rho_bar = -(rho0/g)*N2
 	!
 	!  FOR NOW, AT MOST EXPAND ABOUT A LINEAR PROFILE WHEN USING LAPLACIAN DIFFUSION,
 	!  i.e. SECOND DERIV=0
 	!  **** THE ALGORITHM WILL NOT APPLY DIFFUSION TO THE AMBIENT PROFILES ****
 	!  FOR HIGHER ORDER OPERATORS, MAKE SURE THAT the 2p'th DERIV OF S1_BAR=0
 	!-------------------------------------------------------------------------
 	s1_bar(:,1) = -(rho0/g)*N*N*( z(:) - Lz )   ! s1_bar(z)     [kg/m3]
 	s1_bar(:,2) = -(rho0/g)*N*N                 ! d/dz s1_bar   [kg/m4]
 	s1_bar(:,3) =  0.d0                         ! d2/dz2 s1_bar [kg/m5]
 	ambient_profile(1) = .FALSE.                ! if deriv is nonzero, this needs to be TRUE
 
 	if( do_second_scalar ) then
 		s2_bar(:,1) = 0.d0    ! s2_bar(z)
 		s2_bar(:,2) = 0.d0    ! d/dz s2_bar
 		s2_bar(:,3) = 0.d0    ! d2/dz2 s2_bar
 		ambient_profile(2) = .FALSE.             ! if deriv is nonzero, this needs to be TRUE
 	endif

 return
end subroutine prescribe_ambient_scalar_profiles



subroutine user_ics(x,y,z,F,id,nx,ny,nz)
	!!-----------------------------------------------------------
	!!  Inputs and outputs are all in dimensional units.
	!!
	!!  inputs:
	!!   nx,ny,nz  array size parameters
	!!   x,y,z     coordinate arrays (in meters) that define
	!!             the region over which forcing functions
	!!             are to be evaluated
	!!
	!!   id        field id [u,v,w,s1',s2'] <---[1,2,3,4,5]
	!!             where s1 = s1bar(z) + s1'(y,x,z,t)
	!!             and   s2 = s2bar(z) + s2'(y,x,z,t)
	!!
	!!  outputs:
	!!          F(j,i,k)  i.e. F(y,x,z)
	!!          the initial conditions for variable id
	!!          F is to be returned in DIMENSIONAL units
	!!          (i.e. m/s or deg C etc )
	!!
	!!  N.B. STORAGE INDEX CONVENTION: F(y,x,z)
	!!
	!!-----------------------------------------------------------
 	use mpi_params,             only: myid,comm,ierr,numprocs
 	use decomposition_params,   only: proc_row,proc_col,YBLOCK,p1,p2
 	use independent_variables,  only: Lx,Ly,Lz
 	use dependent_variables,    only: s1,s1_bar
 	use methods_params,         only: restart
 	use user_params,            only: change_default_values,eval_rho_bar,  \
 	                                  delta_U,h_u,h_rho,z0_u,z0_rho,z_off,amp_pert


 	implicit none
 	include 'mpif.h'
 	integer, intent(in)                  ::  id,nx,ny,nz	
 	real(kind=8), intent(in)             ::  x(nx),y(ny),z(nz)  ! local portions
 	real(kind=8), intent(out)            ::  F(ny,nx,nz)
 	integer                              ::  i,j,k

 	!------------------------------------------------------------
 	! user's parameter declarations
 	!------------------------------------------------------------
 	real(kind=8)              :: pi,urv,A,rho_bar
 	real(kind=8), external    :: mytanh
  	
  	
  	! Change any default values specified in user_params_module
  	call change_default_values
  	
  	! If the run is a restart; return.
  	if( restart ) return           

  	pi = 4.d0 * atan(1.d0)
  	z0_u = Lz/2.d0
  	z0_rho = z0_u-z_off


	!------------------------------------------------------------
	! Code section for prescribing initial conditions at t = 0,
	!  Evaluate large scale IW defined in user_params within
	!  the computational domain. XVAL, YVAL and ZVAL are coord
	!  in the larger, outer domain where the IW mode is defined.
	!------------------------------------------------------------
	select case(id)
	
    case(1)        !! ics for u [m/s]
		do k=1,nz
      		do i=1,nx
       			do j=1,ny
       				F(j,i,k) = 0.d0
       			enddo
       		enddo
       	enddo
     
	case(2)        !! ics for v [m/s]
		do k=1,nz
      		do i=1,nx
       			do j=1,ny
       				F(j,i,k) = (delta_U/2.d0)*mytanh( 2.d0*(z(k) - z0_u)/h_u )
       			enddo
       		enddo
       	enddo 

	case(3)        !! ics for w [m/s]
		do k=1,nz
      		do i=1,nx
       			do j=1,ny
       				F(j,i,k) = 0.d0
       			enddo
       		enddo
       	enddo

    case(4)        !! ics for s1' = density [kg/m3] with s1_bar=0
    	do k=1,nz
    		rho_bar = eval_rho_bar(z(k))
      		do i=1,nx
       			do j=1,ny
       			
       				call RANDOM_NUMBER( urv )                 ! in [0,1]
    				urv = 2.d0*(urv-0.5d0)                    ! in [-1 1]
    				A = amp_pert*(1.d0 + 1.d-3*urv)           ! [kg/m3] 
    				
       				F(j,i,k) =  rho_bar + A*(cos(2.d0*pi*y(j)/Ly))**2 / ( cosh( 2.d0*(z(k) - z0_rho)/h_rho ) )**2
       				
       			enddo
       		enddo
       	enddo
  

    case(5)        !! ics for s2'
		F(:,:,:) = 0.d0

    case default
      stop 'error in calling user_ics: id out of range'
	end select
	call mpi_barrier(comm,ierr)
 return
end subroutine user_ics


subroutine user_forcing(x,y,z,t,F,id,nx,ny,nz,call_again)
!-------------------------------------------------------------------                        
!  User defined subroutine to add time/space dependent
!  forcing to the equations of motion. 
!
!
!  inputs:
!   nx,ny,nz  array size parameters
!   x,y,z     coordinate arrays [m] that define
!             the grid points at which the forcing
!             functions are to be evaluated
!   t         current time in seconds
!   id        field id [u,v,w,s1',s2'] <---[1,2,3,4,5]
!
!  outputs:
!            F(j,i,k) i.e. F(y,x,z)
!            the rhs term for equation id
!            F is to be returned in DIMENSIONAL units
!            (i.e. m/s2 or deg C/s etc )
!
!            call_next_time
!            set false if field id is not being forced
!            e.g. suppose do_forcing=.true., but only u
!            is actually forced, set call_next_time=.FALSE.
!            for v,w,s1 & s2 to prevent adding zeros to
!            the respective rhs's at each time step.
!
!  NB STORAGE INDEX CONVENTION: F(y,x,z)
!
!-------------------------------------------------------------------
	use mpi_params,             only: myid
	use dimensional_scales,     only: g,rho0,f0,nu,kappa
	use independent_variables,  only: Lx,Ly,Lz,dt,t0,tf 
	use dependent_variables,    only: u,v,w,s1,s2,s1_bar,s2_bar
	use decomposition_params,   only: proc_row,proc_col,YBLOCK,p1,p2
	implicit none
	integer                                :: i,j,k
	integer                                :: id,nx,ny,nz
	real(kind=8)                           :: x(nx),y(ny),z(nz),t
	real(kind=8)                           :: F(ny,nx,nz)
	logical                                :: call_again
	real(kind=8),save                      :: pi,z0_u,z0_rho
	logical,save                           :: first_entry=.TRUE.
	include 'mpif.h'

	!------------------------------------------------------------
	! user's parameter declarations
	!------------------------------------------------------------
	real(kind=8)                          :: dx=1.d0, dy=1.d0, dz=1.d0

	!------------------------------------------------------------
	! first entry calculation of saved parameters
	!------------------------------------------------------------
	if(first_entry) then
 		pi=4.d0*atan(1.d0)    
		if(nx > 1) dx = x(2)-x(1)                     ! [m]
		if(ny > 1) dy = y(2)-y(1)                     ! [m]
		if(nz > 1) dz = z(2)-z(1)                     ! [m]   
						
		!---------------------------------------------------------------------
		!  write a simple text file with some parameter values...
		!---------------------------------------------------------------------  
		if( myid==0 ) then
			open(1,file='output/simulation_params')
			write(1,*) 'p1 x p2                ', p1,p2
			write(1,*) 'tf - t0 [hours]        ',(tf-t0)/3600.d0
			write(1,*) ' '
			write(1,*) 'Lx [m]                 ', Lx
			write(1,*) 'Ly [m]                 ', Ly
			write(1,*) 'Lz [m]                 ', Lz
			write(1,*) 'dx [m]                 ', dx
			write(1,*) 'dy [m]                 ', dy
			write(1,*) 'dz [m]                 ', dz
			write(1,*) ' '
			write(1,*) 'dt [s]                 ', dt
			write(1,*) 'nu [m/s2]              ', nu
			write(1,*) 'kappa [m/s2]           ', kappa(1)
			write(1,*) 'dx^2/nu*dt             ', dx*dx/(nu*dt)
			write(1,*) ' '
			close(1)
		endif
		     
		first_entry=.FALSE.
	endif 
   

	if( id == 1 ) then       !  F <===  specify additional forcing on rhs of u equation 
		do k=1,nz
			do i=1,nx
				do j=1,ny
     				F(j,i,k) = 0.d0
    			enddo
   			enddo
  		enddo  
  		call_again = .FALSE.  
  		
	elseif( id==2 ) then     !  F <===  specify additional forcing on rhs of v equation
		do k=1,nz
			do i=1,nx
				do j=1,ny
     				F(j,i,k) =  0.d0
    			enddo
   			enddo
  		enddo  
  		call_again = .FALSE.
	
	elseif( id==3 ) then     !  F <===  specify additional forcing on rhs of w equation
		do k=1,nz
			do i=1,nx
				do j=1,ny
     				F(j,i,k) = 0.d0  
    			enddo
   			enddo
  		enddo  
  		call_again = .FALSE.
	
	elseif( id==4 ) then      !  F <===  specify additional forcing on rhs of s1 equation
		do k=1,nz			
			do i=1,nx
				do j=1,ny
					F(j,i,k) = 0.d0	
    			enddo
   			enddo
  		enddo  
  		call_again = .FALSE.
	
	elseif( id==4 ) then      !  F <===  specify additional forcing on rhs of s2 equation
		do k=1,nz
			do i=1,nx
				do j=1,ny
     				F(j,i,k) = 0.d0  ! or whatever you want to add
    			enddo
   			enddo
  		enddo  
  		call_again = .FALSE.   
	endif
	
 return  
end subroutine user_forcing

	
subroutine user_bvals_ew    
	!-------------------------------------
	! routine must be here, even if empty
	!-------------------------------------
	implicit none
	return
end subroutine user_bvals_ew
	
subroutine user_bvals_ns    
	!-------------------------------------
	! routine must be here, even if empty
	!-------------------------------------
	implicit none
	return
end subroutine user_bvals_ns
	
subroutine user_bvals_bt    
	!-------------------------------------
	! routine must be here, even if empty
	!-------------------------------------
	implicit none
	return
end subroutine user_bvals_bt
