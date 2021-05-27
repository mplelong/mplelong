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
 

 !-------------------------------------------------------------------------
 ! For this problem, expand around linear density profile with specified N
 !  N2 = -(g/rho0)* d/dz rho_bar  ==> d/dz rho_bar = -(rho0/g)*N2
 !
 !  FOR NOW, AT MOST EXPAND ABOUT A LINEAR PROFILE WHEN USING LAPLACIAN DIFFUSION,
 !  i.e. SECOND DERIV=0
 !  **** THE ALGORITHM WILL NOT APPLY DIFFUSION TO THE AMBIENT PROFILES ****
 !  FOR HIGHER ORDER OPERATORS, MAKE SURE THAT the 2p'th DERIV OF S1_BAR=0
 !-------------------------------------------------------------------------
 s1_bar(:,1) = 0.d0                          ! s1_bar(z)     [kg/m3]
 s1_bar(:,2) = 0.d0                          ! d/dz s1_bar   [kg/m4]
 s1_bar(:,3) = 0.d0                          ! d2/dz2 s1_bar [kg/m5]
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
 	use user_params,            only: change_default_values,rho_north,rho_south,   &
 	                                  delta_rho_north,gamma_ns,xnoise,pi,NS_profiles

 	implicit none
 	include 'mpif.h'
 	integer, intent(in)                  ::  id,nx,ny,nz	
 	real(kind=8), intent(in)             ::  x(nx),y(ny),z(nz)  ! local portions
 	real(kind=8), intent(out)            ::  F(ny,nx,nz)
 	integer                              ::  i,j,k
 	logical, save                        :: first_entry = .TRUE.

 	!------------------------------------------------------------
 	! user's parameter declarations
 	!------------------------------------------------------------
 	real(kind=8)              :: urv
 	real(kind=8)              :: y0,gamma,xx           

	
	! If the run is a restart; return.
  	if( restart ) return


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
       				F(j,i,k) = 0.d0
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
    	y0 = Ly/2.d0
    	gamma = Ly*gamma_ns
    	do k=1,nz
      		do i=1,nx
       			do j=1,ny
       			
       				F(j,i,k)  = rho_south(k)*( .5*(1.-tanh((y(j)-y0)/gamma)) )    &
                   	          + rho_north(k)*( .5*(1.+tanh((y(j)-y0)/gamma)) )
                                               
        			! add just a bit of noise to seed instabilities
        			call RANDOM_NUMBER( urv )                 ! in [0,1]
        			urv = 2.d0*(urv-0.5d0)                    ! in [-1 1] 
        
        			! confine noise to interior of the domain
        			xx = sin( pi*z(k)/Lz )**2         
       				xx = xx*sin( pi*x(i)/Lx )**2
        			xx = xx*sin( pi*y(j)/Ly )**2
        			F(j,i,k) = F(j,i,k) + xnoise*delta_rho_north*urv*xx
        
       			enddo
       		enddo
       	enddo
  

    case(5)        !! ics for s2'
		do k=1,nz
      		do i=1,nx
       			do j=1,ny
       				F(j,i,k) = 0.d0
       			enddo
       		enddo
       	enddo

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
	use mpi_params,             only: myid,comm,ierr
	use dimensional_scales,     only: g,rho0,f0
	use independent_variables,  only: Lx,Ly,Lz,dt,t0,tf 
	use dependent_variables,    only: u,v,w,s1,s2,s1_bar,s2_bar
	use decomposition_params,   only: proc_row,proc_col,YBLOCK,p1,p2	
	use user_params
	
	implicit none
	integer                                :: i,j,k
	integer                                :: id,nx,ny,nz
	real(kind=8)                           :: x(nx),y(ny),z(nz),t
	real(kind=8)                           :: F(ny,nx,nz)
	logical                                :: call_again
	real(kind=8)                           :: urv
	real(kind=8), external                 :: myexp
	logical,save                           :: first_entry=.TRUE.
	include 'mpif.h'

	!----------------------------------------------------------------------------
	! user's parameter declarations
	!----------------------------------------------------------------------------
	real(kind=8)                           :: dx=1.d0, dy=1.d0, dz=1.d0
	real(kind=8), save                     :: beta,f_tilde,y0,xx
	character(len=80)                      :: cwd
	character*(MPI_MAX_PROCESSOR_NAME)     :: hostname
	integer                                :: len
	real(kind=8)                           :: tau,taux,tauy
	real(kind=8),save                      :: cd_wind,rho_air
	real(kind=8)                           :: u10_of_t,theta_of_t
	real(kind=8)                           :: t_left,t_right
	integer                                :: i0,i1

	!------------------------------------------------------------
	! first entry calculation of saved parameters
	!------------------------------------------------------------
	if(first_entry) then
	
		call getcwd(cwd)    ! gnu intrinsic.... ok on other compilers? intel ok
		call mpi_get_processor_name(hostname,len,ierr)
		
		if(nx > 1) dx = x(2)-x(1)                 	! [m]
		if(ny > 1) dy = y(2)-y(1)					! [m]
		if(nz > 1) dz = z(2)-z(1)					! [m] 
		
		y0 = Ly/2.d0                                ! [m]
		beta = 2*Omega*cos(lat)/R_earth				! [1/ms]
  		f_tilde = 2*Omega*cos(lat)					! [1/s]
  		if( .NOT. beta_plane ) beta = 0.d0
  		if( .NOT. NT_terms   ) f_tilde = 0.d0 
  		
  		cd_wind = 1.2d-3                   ! [1]
		rho_air = 1.22                     ! [kg/m3] 
				
		allocate( y_relax_window(ny,2) )  
  		do j=1,ny
   			y_relax_window(j,1) = myexp( -((y(j)-0.d0)/(sponge_scale*1000.d0))**2 )   ! 1 near y=0
   			y_relax_window(j,2) = myexp( -((y(j)- Ly )/(sponge_scale*1000.d0))**2 )   ! 1 near y=Ly
  		enddo
  		
  		allocate( drag_window(nz) )
		do k=1,nz
			drag_window(k) = myexp(-( z(k)/delta_bottom_drag)**2)
		enddo
		if( .NOT. bottom_drag   ) drag_window(:) = 0.d0
		
		allocate ( tau_window(nz) )    ! assume we keep z_MLB constant, easy to generalize
		do k=1,nz
			tau_window(k) = myexp( -( (z(k)-Lz)/depth_MLB_initial)**2 )  ! [1]  confined to surface mixed layer
		enddo
		
		if( wind_stress ) then
			call compute_wind_speed_direction   ! compute or read and store random sequences
		else
			tau_window(:) = 0.d0
		endif
		
		
		!----------------------------------------------------------------------------------------
		!  write a simple text file with some parameter values...
		!----------------------------------------------------------------------------------------  
		if( myid==0 ) then
			open(1,file='output/simulation_params')
			write(1,*) 'hostname                          ',trim(hostname)
			write(1,*) 'flow_solve root dir               ',trim(cwd)
			write(1,*) 'p1 x p2                           ', p1,p2
		    write(1,*) 'tf - t0 [hours]                   ',(tf-t0)/3600.
			write(1,*) ' '
			write(1,*) 'Lx [km]                           ', Lx/1000.
			write(1,*) 'Ly [km]                           ', Ly/1000.
			write(1,*) 'Lz [m]                            ', Lz
			write(1,*) 'dx [m]                            ', dx
			write(1,*) 'dy [m]                            ', dy
			write(1,*) 'dz [m]                            ', dz
			write(1,*) ' '
			write(1,*) 'dt [s]                            ', dt
			write(1,*) 'latitude [deg]                    ', lat*360./(2.*pi)
			write(1,*) 'dt wind data [s]                  ', dt_rand
			write(1,*) 'f0 [1/s]                          ', f0
			write(1,*) 'f~ [1/s]                          ', f_tilde
			write(1,*) 'beta [1/s]                        ', beta
			write(1,*) '2pi/f0 [hours]                    ', (2.*pi/f0)/3600.
			write(1,*) ' '
			write(1,*) 'Non-traditional                   ', NT_terms
			write(1,*) 'beta plane                        ', beta_plane
			write(1,*) 'bottom drag                       ', bottom_drag
			write(1,*) 'wind stress                       ', wind_stress
			write(1,*) 'read wind data?                   ', read_wind_data 
			write(1,*) ' '     
			write(1,*) 'u10 mean [m/s]                    ', u10_mean
			write(1,*) 'u10 st dev [m/s]                  ', u10_stdev
			write(1,*) 'u10 relax time [hrs]              ', u10_tau/3600.
			write(1,*) ' '
			write(1,*) 'NS surface delta rho [kg/m3]      ', delta_rho_surf
			write(1,*) 'top to bottom delta rho N [kg/m3] ', delta_rho_north
			write(1,*) 'rho decay scale [m]               ', lambda_north
			write(1,*) 'depth scale NS diff [m]           ',lambda_star 
			write(1,*) 'N/S relaxation time scale [hrs]   ',t_damp /3600.	
			close(1)
		endif
      
		call mpi_barrier(comm,ierr)
		first_entry=.FALSE.
	endif 
	
	
	
	
	!---------------------------------------------------------------------------
 	! get the x and y components of the stochastic wind stress at current time 
 	!---------------------------------------------------------------------------
	if( wind_stress .and. id < 3 ) then
		!---------------------------------------------------------------------------
 		!   interpolate between stored,random values to get 
		!   wind speed and direction at the current time
		!---------------------------------------------------------------------------
		i0 = floor( (t-t0) /dt_rand ) + 1      ! index of coarse time just lt t
		i1 = i0 + 1                            ! index of coarse time just gt t
		t_left = t0 + (i0-1)*dt_rand           ! coarse time in [s] just lt t
		t_right = t_left + dt_rand             ! coarse time in [s] just gt t
 
		u10_of_t = u10(i0) + ((t-t_left)/dt_rand)*( u10(i1)-u10(i0) )
		theta_of_t = theta(i0) + ((t-t_left)/dt_rand)*( theta(i1)-theta(i0) )
 
 		if( t < t_left .OR. t > t_right ) then
 			write(0,*) 'Error aligning wind data time grid: t,dt_rand,i0,i1 ',t,dt_rand,i0,i1
  			STOP
 		endif
 		
 		!---------------------------------------------------------------------------
 		!   convert to x and y components of surface wind stress 
 		!---------------------------------------------------------------------------
 		tau = cd_wind * rho_air * u10_of_t**2     	! [kg/ms2]
 		taux = tau * cos( theta_of_t )				! x component
 		tauy = tau * sin( theta_of_t )				! y component	
 	else 	
 		taux = 0.d0
 		tauy = 0.d0	
 	endif

	
	   

	if( id == 1 ) then       !  F <===  specify additional forcing on rhs of u equation 
		do k=1,nz
			do i=1,nx
				do j=1,ny
				
					!-----------------------------------------------------
					! beta effect + NT term
					!-----------------------------------------------------
     				F(j,i,k) = beta*(y(j)-y0)*v(j,i,k) - f_tilde*w(j,i,k) 	
     				
     				!--------------------------------------------------------------- 
     				! Rayleigh damp toward zero near N and S edges of domain 
     				!--------------------------------------------------------------- 
     				F(j,i,k) = F(j,i,k)                                      &
        				- y_relax_window(j,1)*(1./t_damp)*(u(j,i,k) - 0.d0)  &
        				- y_relax_window(j,2)*(1./t_damp)*(u(j,i,k) - 0.d0)
        				
        			!---------------------------------------------------------------------- 
     				! quadratic bottom drag    
     				!----------------------------------------------------------------------     
					xx = CD * sqrt( u(j,i,k)**2 + v(j,i,k)**2 )/delta_bottom_drag  ! [1/s]     
					F(j,i,k) = F(j,i,k) - drag_window(k)*xx*u(j,i,k)               ! [m/s2]
					
					!---------------------------------------------------------------------- 
     				! wind stress
     				!----------------------------------------------------------------------
     				xx = (1.d0/rho0) * taux/depth_MLB_initial                	   ! [m/s2] 
     				F(j,i,k) = F(j,i,k) + tau_window(k)*xx
        				
    			enddo
   			enddo
  		enddo   
  		call_again = .TRUE.  
  		
	elseif( id==2 ) then     !  F <===  specify additional forcing on rhs of v equation
		do k=1,nz
			do i=1,nx
				do j=1,ny
				
					!-----------------------------------------------------
					! beta effect 
					!-----------------------------------------------------
     				F(j,i,k) = -beta*(y(j)-y0)*u(j,i,k)
     				
     				!--------------------------------------------------------------- 
     				! Rayleigh damp toward zero near N and S edges of domain 
     				!--------------------------------------------------------------- 
     				F(j,i,k) = F(j,i,k)                                      &
        				- y_relax_window(j,1)*(1./t_damp)*(v(j,i,k) - 0.d0)  &
        				- y_relax_window(j,2)*(1./t_damp)*(v(j,i,k) - 0.d0)
        				
        			!---------------------------------------------------------------------- 
     				! quadratic bottom drag    
     				!----------------------------------------------------------------------     
					xx = CD * sqrt( u(j,i,k)**2 + v(j,i,k)**2 )/delta_bottom_drag  ! [1/s]     
					F(j,i,k) = F(j,i,k) - drag_window(k)*xx*v(j,i,k)               ! [m/s2]
					
					!---------------------------------------------------------------------- 
     				! wind stress
     				!----------------------------------------------------------------------
     				xx = (1.d0/rho0) * tauy/depth_MLB_initial                	   ! [m/s2] 
     				F(j,i,k) = F(j,i,k) + tau_window(k)*xx
        										
    			enddo
   			enddo
  		enddo  
  		call_again = .TRUE.
	
	elseif( id==3 ) then     !  F <===  specify additional forcing on rhs of w equation
		do k=1,nz
			do i=1,nx
				do j=1,ny
				
					!-----------------------------------------------------
					! NT term 
					!-----------------------------------------------------
     				F(j,i,k) = f_tilde * u(j,i,k)
     				
     				!---------------------------------------------------------------
     				! Rayleigh damp toward zero near N and S edges of domain 
     				!--------------------------------------------------------------- 
     				F(j,i,k) = F(j,i,k)                                      &
        				- y_relax_window(j,1)*(1./t_damp)*(w(j,i,k) - 0.d0)  &
        				- y_relax_window(j,2)*(1./t_damp)*(w(j,i,k) - 0.d0)
        						
    			enddo
   			enddo
  		enddo  
  		call_again = .TRUE.
	
	elseif( id==4 ) then      !  F <===  specify additional forcing on rhs of s1 equation
		do k=1,nz
			do i=1,nx
				do j=1,ny
     				
     				!---------------------------------------------------------------------- 
     				! Rayleigh damp toward specified profiles near N and S edges of domain 
     				!---------------------------------------------------------------------- 
     				F(j,i,k) = - y_relax_window(j,1)*(1./t_damp)*(s1(j,i,k) - rho_south(k))   &
     						   - y_relax_window(j,2)*(1./t_damp)*(s1(j,i,k) - rho_north(k))       
       
    			enddo
   			enddo
  		enddo  
  		call_again = .TRUE.
	
	elseif( id==5 ) then      !  F <===  specify additional forcing on rhs of s2 equation
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

subroutine user_bvals_EW(x,y,z,t,id,VALS,DERIVS,nx,ny,nz,side,Lx)
	!--------------------------------------------------------------------------
	! User routine to supply externally obtained boundary values and normal
	! derivatives at the EAST and WEST boundaries x=0 and x=Lx at
	! the y and z values requested at time t. All values in dimensional units.
	! "side" is either 'E' for x=0 or 'W' for x=Lx requested values.
	!--------------------------------------------------------------------------
	
	!-----------------------------------------------------------------------------
	!   These are user defined functions that return the actual boundary values.
	!   Generally, the user will read and store boundary values coarser in
	!   time than needed, and then interpolate to the desired time, returning the
	!   requested values via a user defined function. Here we just evaluate an analytical
	!   external solution directly at the locations and times requested by the
	!   calling routine.
	!-----------------------------------------------------------------------------
	use user_params
	!-----------------------------------------------------------------------------
		
	implicit none
	integer, intent(in)                    :: id                  ! 1,2,3,4,5 for u,v,w,s1,s2
	integer, intent(in)                    :: nx,ny,nz            ! size of local portions of domain
	real(kind=8), intent(in)               :: x(nx),y(ny),z(nz)   ! local values of y and z
	real(kind=8), intent(out)              :: VALS(ny,nz)
	real(kind=8), intent(out)              :: DERIVS(ny,nz)
	real(kind=8), intent(in)               :: t,Lx                ! time in [s], global size Lx in [m]
	character(len=80),intent(in)           :: side 
		
 return	          
end subroutine user_bvals_EW

subroutine user_bvals_NS(x,y,z,t,id,VALS,DERIVS,nx,ny,nz,side,Ly)
	!--------------------------------------------------------------------------
	! User routine to supply externally obtained boundary values and normal
	! derivatives at the NORTH and SOUTH boundaries y=0 and y=Ly at
	! the x and z values requested at time t. All values in dimensional units.
	! "side" is either 'S' for y=0 or 'N' for y=Ly requested values.
	!--------------------------------------------------------------------------
	
	!-----------------------------------------------------------------------------
	!   These are user defined functions that return the actual boundary values.
	!   Generally, the user will read and store boundary values coarser in
	!   time than needed, and then interpolate to the desired time, returning the
	!   requested values via a user defined function. Here we just evaluate an analytical
	!   external solution directly at the locations and times requested by the
	!   calling routine.
	!-----------------------------------------------------------------------------
	use user_params,                    only: rho_north,rho_south
	!-----------------------------------------------------------------------------
		
	implicit none
	integer, intent(in)                    :: id                  ! 1,2,3,4,5 for u,v,w,s1,s2
	integer, intent(in)                    :: nx,ny,nz            ! size of local portions of domain
	real(kind=8), intent(in)               :: x(nx),y(ny),z(nz)   ! local values of y and z
	real(kind=8), intent(out)              :: VALS(nx,nz)
	real(kind=8), intent(out)              :: DERIVS(nx,nz)
	real(kind=8), intent(in)               :: t,Ly                ! time in [s], global size Ly in [m]
	character(len=80),intent(in)           :: side 
	character(len=80), save                :: dir='y'             ! return y derivs at N/S bdries 
	integer                                :: i,k
	real(kind=8)                           :: xval,yval,zval
	logical                                :: do_south,do_north   

	do_south=.FALSE. ;  do_north=.FALSE.
	if(side=='S' .and. y(1) == 0.d0 )   do_south = .TRUE.
	if(side=='N' .and. y(ny) == Ly )	do_north = .TRUE.
		
	if( do_south )  then
		if( id==4 ) then
			do k=1,nz
		 		VALS(:,k) = rho_south(k)  ! set bval to prescribed south profile
			enddo
		else
			VALS(:,:) = 0.d0              ! all other variables = 0 at outer edge of sponge
		endif
		DERIVS(:,:) = 0.d0                ! all derivs d/dy=0 at N/S boundaries
	endif
	
	if( do_north )  then
		if( id==4 ) then
			do k=1,nz
		 		VALS(:,k) = rho_north(k)  ! set bval to prescribed north profile
			enddo
		else
			VALS(:,:) = 0.d0              ! all other variables = 0 at outer edge of sponge
		endif
		DERIVS(:,:) = 0.d0                ! all derivs d/dy=0 at N/S boundaries
	endif
 return	          
end subroutine user_bvals_NS

subroutine user_bvals_BT(x,y,z,t,id,VALS,DERIVS,nx,ny,nz,side,Lz)
	!--------------------------------------------------------------------------
	! User routine to supply externally obtained boundary values and normal
	! derivatives at the BOTTOM and TOP boundaries z=0 and z=Lz at
	! the x and y values requested at time t. All values in dimensional units.
	! "side" is either 'B' for z=0 or 'T' for z=Lz requested values.
	!--------------------------------------------------------------------------
	
	!-----------------------------------------------------------------------------
	!   These are user defined functions that return the actual boundary values.
	!   Generally, the user will read and store boundary values coarser in
	!   time than needed, and then interpolate to the desired time, returning the
	!   requested values via a user defined function. Here we just evaluate an analytical
	!   external solution directly at the locations and times requested by the
	!   calling routine.
	!-----------------------------------------------------------------------------
	use user_params
	!-----------------------------------------------------------------------------
		
	implicit none
	integer, intent(in)                    :: id                  ! 1,2,3,4,5 for u,v,w,s1,s2
	integer, intent(in)                    :: nx,ny,nz            ! size of local portions of domain
	real(kind=8), intent(in)               :: x(nx),y(ny),z(nz)   ! local values of y and z
	real(kind=8), intent(out)              :: VALS(nx,ny)
	real(kind=8), intent(out)              :: DERIVS(nx,ny)
	real(kind=8), intent(in)               :: t,Lz                ! time in [s], global size Ly in [m]
	character(len=80),intent(in)           :: side 
	
 return	          
end subroutine user_bvals_BT

subroutine compute_wind_speed_direction
	!----------------------------------------------------------------
	!  Routine to compute the random 10 m windspeed and direction
	!  to be used in the surface forcing
	!  input 
	!       x,y,z  [m]  1d arrays of length nx,ny,nz
	!       t      [s]  current time
	!       dt     [1]  D'LESS time step
	!       myid   processor id
	!  output 
	!  z_MLB(y,x) for processor myid in [m]  top plane of YBLOCK format
	!------------------------------------------------------------------
	use user_params,            only: u10,theta,u10_mean,u10_stdev,u10_tau,  &
									  dt_rand,read_wind_data,t_offset
	use mpi_params,             only: myid,comm,ierr
	use independent_variables,  only: dt,t0,tf         ! dt, t0,tf in [s]

	implicit none
	real(kind=8), external         :: norm_rand
	real(kind=8)                   :: pi,mean,std_dev,t,ws,wd
	real(kind=8)                   :: mu,sigma,tau,xx
	real(kind=8),allocatable       :: tmp(:)
	integer                        :: i,line,root,count,n_seq
	character*80                   :: method='Ornstein-Uhlenbeck'
	include 'mpif.h'
                            
	pi = 4.d0*atan(1.d0)
  
  
	n_seq = (tf-t0) / dt_rand + 5                         ! extend outside range for interpolation
	allocate( u10(n_seq),theta(n_seq), tmp(n_seq) )
  
  
	if( myid==0 .AND. .NOT. read_wind_data ) then         ! create random sequences
  
   
		if( method == 'normal_dist' ) then
			!--------------------------------------------------------------  
			!  wind speed normally distributed with specified mean,variance  
			!--------------------------------------------------------------
			mean = u10_mean
			std_dev = u10_stdev 
			do i=1,n_seq
				u10(i) =  norm_rand(mean,std_dev)   ! gaussian with specified mean and standard deviation
			enddo
   	
			!------------------------------------------------------  
			!  wind direction uniformly distributed in -pi to pi  
			!------------------------------------------------------
			do i=1,n_seq
				call random_number(theta(i))             !  theta in [0,1]
				theta(i) = 2.d0*(theta(i) - 0.5d0)       !  theta in [-1,1]
				theta(i) = pi*theta(i)                   !  theta in [-pi,pi]
			enddo
   	
		elseif( trim(method) == 'Ornstein-Uhlenbeck' ) then
			!--------------------------------------------------------------  
			!  wind speed given by Ornstein-Uhlenbeck process w/
			!  mean mu=u10_mean and st. dev. sigma=u10_stdev 
			!  and tau = u10_tau
			!--------------------------------------------------------------
			mean = 0.d0
			std_dev = 1.d0 
			mu = u10_mean
			sigma = u10_stdev
			tau = u10_tau
			xx=sigma*sqrt(2.d0/tau)*sqrt(dt_rand)
			u10(1) = mu  !mean
			do i=2,n_seq
				u10(i) = u10(i-1) + dt_rand*( -(u10(i-1)-mu)/tau ) + xx*norm_rand(mean,std_dev)
			enddo
    
			!--------------------------------------------------------------  
			!  wind direction given by Ornstein-Uhlenbeck process w/
			!  mean mu=pi/4. and st. dev. sigma=pi/8
			!  and tau = u10_tau 
			!--------------------------------------------------------------
			mean = 0.d0
			std_dev = 1.d0 
			mu = pi/4.d0
			sigma = pi/8.d0
			tau = u10_tau
			xx=sigma*sqrt(2.d0/tau)*sqrt(dt_rand)
			theta(1) = mu  ! mean
			do i=2,n_seq
				theta(i) = theta(i-1) + dt_rand*( -(theta(i-1)-mu)/tau ) + xx*norm_rand(mean,std_dev)
			enddo
    
		endif
   
     
	elseif( myid==0 .AND. read_wind_data ) then       ! read sequences from file
  
		i=0
		t = -999999.999    ! just to start the loop
   
		open(1,file='input/wind_speed_direction')
   
		do while ( t <= tf + t_offset )      
			read(1,*) t, ws, wd       !  [s],[m//s],[radians]            
			if( t >= t0 + t_offset ) then
				i = i + 1
				u10(i) = ws
				theta(i) = wd
			endif    
		enddo
      
	   close(1)
    
	endif
  
  
	!--------------------------------------------------------------------------------------
	! keep and use the myid=0 sequence for each processor, i.e. discard others
	!--------------------------------------------------------------------------------------
	root = 0
	count = n_seq
  
	call MPI_BARRIER(comm,ierr)
	call MPI_Bcast(u10,count,MPI_DOUBLE_PRECISION,root,comm,ierr)
	call MPI_Bcast(theta,count,MPI_DOUBLE_PRECISION,root,comm,ierr)

  
	!--------------------------------------------------------------------------------------
	! write out the wind data to a text file
	!--------------------------------------------------------------------------------------  
	if( myid==0 ) then
		open(1,file='output/1D/wind_speed_direction')
			do i=1,n_seq
				write(1,200) (i-1.d0)*dt_rand,u10(i),theta(i)
			enddo
		close(1)
	endif
200 FORMAT (3e15.6)  
 
	deallocate( tmp )
	call MPI_BARRIER(comm,ierr)
 return
end subroutine compute_wind_speed_direction




real(kind=8) function norm_rand(mean, std_dev)

	!The Marsalgia Polar Method
	!This method (or algorithm) is a more efficient version of the more commonly known Box-Muller transform. 
	!The Marsalgia method is more efficient because it uses fewer calls to trigonometric functions, which 
	!are computationally expensive to run, replacing these calls with a rejection procedure. 
	! uses the marsaglia polar method to create pseudo-random number pairs that
	! have a normal distribution. The function saves one of the random numbers from
	! the generated pair as a spare to be returned for the next call to the function.


	! Returns a real(kind=8) scalar
	! uses f90 intrinsic routine random_number

	!  KBW:  made all reals (kind=8), all constants double precision
    implicit none
    real(kind=8), intent(in) :: mean, std_dev
    real(kind=8) :: x, y, r
    real(kind=8), save :: spare
    logical, save :: has_spare
    
    ! use a spare saved from a previous run if one exists
    if (has_spare) then
        has_spare = .FALSE.
        norm_rand = mean + (std_dev * spare)
        return
    else
        r = 1.d0
        do while ( r >= 1.d0 )
            ! generate random number pair between 0 and 1
            call random_number(x)
            call random_number(y)
            ! normalise random numbers to be in square of side-length = R
            x = (x * 2.d0) - 1.d0
            y = (y * 2.d0) - 1.d0
            r = x*x + y*y
        end do
        
        ! calculate the co-efficient to multiply random numbers x and y
        ! by to achieve normal distribution
        r = sqrt((-2.d0 * log(r)) / r)
        norm_rand = mean + (std_dev * x * r)
        spare = y * r
        has_spare = .TRUE.
        return
    end if
end function norm_rand


