subroutine explicit_step
!-------------------------------------------------------------------------
!  Routine to take an explicit time step via Euler,
!  2nd or third order Adams Bashforth method. The differential
!  equations are as follows:
!
!  d/dt u   = rhs1
!  d/dt v   = rhs2
!  d/dt w   = rhs3
!  d/dt s1  = rhs4
!  d/dt s2  = rhs5
!
!  notation:
!  dt   time increment
!  MM0   integer index identifying the time t    storage locations
!  MM1   integer index identifying the time t-dt storage locations
!  MM2   integer index identifying the time t-2*dt storage locations
!  MM3   integer index identifying the time t-3*dt storage locations
!
!  explicit_rhs(1,MM0),
!  explicit_rhs(1,MM1),
!  explicit_rhs(1,MM2),
!  explicit_rhs(1,MM3)  
!       ====>   rhs for u at t,t-dt,t-2*dt,,t-3*dt
!  etc for explicit_rhs(2,:),
!          explicit_rhs(3,:),
!          explicit_rhs(4,:) 
!          explicit_rhs(5,:)

!
!  outputs:
!  calculated values of u at time t+dt
!  calculated values of v at time t+dt
!  calculated values of w at time t+dt
!  calculated values of scalar(1) at time t+dt (overwrite)
!  calculated values of scalar(2) at time t+dt (overwrite)
!-------------------------------------------------------------------------
	use mpi_params,              only: myid
	use methods_params,          only: AB_order,do_second_scalar,ambient_profile,high_order_operators
	use decomposition_params
	use independent_variables,   only: dt
	use intermediate_variables,  only: explicit_rhs,ustar,vstar,wstar
	use dimensional_scales,      only: kappa
	use dependent_variables,     only: u,v,w,s1,s2,s1_bar,s2_bar
	use etc
	 
	implicit none
	integer                :: k,k_global
	real(kind=8), save     :: dt2,dt12,dt24
	logical, save          :: first_entry = .TRUE.

	if( first_entry ) then
		dt2 = dt/2.d0
		dt12= dt/12.d0
		dt24= dt/24.d0
		first_entry = .FALSE.
	endif
  
  
	if( step_flag == 'euler') then
		if(myid==0 .and. istep==0) write(0,*) '... 1st order Euler step '
		ustar(:,:,:) = u(:,:,:)  +  dt*explicit_rhs(:,:,:,1,MM0)
		vstar(:,:,:) = v(:,:,:)  +  dt*explicit_rhs(:,:,:,2,MM0)
		wstar(:,:,:) = w(:,:,:)  +  dt*explicit_rhs(:,:,:,3,MM0)
		s1(:,:,:) = s1(:,:,:) + dt*explicit_rhs(:,:,:,4,MM0)
		if( do_second_scalar ) then
			s2(:,:,:) = s2(:,:,:) + dt*explicit_rhs(:,:,:,5,MM0)
		endif
		if(AB_order>1) step_flag='AB2'
  
	elseif( step_flag == 'AB2' ) then
		if(myid==0 .and. istep==1 ) write(0,*) '... 2nd order AB step '

		ustar(:,:,:)  =  u(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,1,MM0)    &
                               -        explicit_rhs(:,:,:,1,MM1) )
                     
		vstar(:,:,:)  =  v(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,2,MM0)    &
                               -        explicit_rhs(:,:,:,2,MM1) )
                     
		wstar(:,:,:)  =  w(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,3,MM0)    &
                               -        explicit_rhs(:,:,:,3,MM1) )
                     
		s1(:,:,:) = s1(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,4,MM0)    &
                               -        explicit_rhs(:,:,:,4,MM1) )
  
		if( do_second_scalar ) then
			s2(:,:,:) = s2(:,:,:)  +  dt2*( 3.d0*explicit_rhs(:,:,:,5,MM0)    &
                                   -        explicit_rhs(:,:,:,5,MM1) )
		endif
		if(AB_order>2) step_flag='AB3'

	elseif( step_flag == 'AB3' ) then
		if(myid==0 .and. istep==2) write(0,*) '... 3rd order AB step '
  
		ustar(:,:,:)  =  u(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,1,MM0)    &
                               - 16.d0*explicit_rhs(:,:,:,1,MM1)            &
                               + 5.d0*explicit_rhs(:,:,:,1,MM2) )
   
  
		vstar(:,:,:)  =  v(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,2,MM0)    &
                               - 16.d0*explicit_rhs(:,:,:,2,MM1)            &
                        	   + 5.d0*explicit_rhs(:,:,:,2,MM2) )
 
		wstar(:,:,:)  =  w(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,3,MM0)    &
                        	   - 16.d0*explicit_rhs(:,:,:,3,MM1)            &
                               + 5.d0*explicit_rhs(:,:,:,3,MM2) )
                    
		s1(:,:,:)  =  s1(:,:,:)  +  dt12*( 23.d0*explicit_rhs(:,:,:,4,MM0)    &
                                 - 16.d0*explicit_rhs(:,:,:,4,MM1)            &
                                 + 5.d0*explicit_rhs(:,:,:,4,MM2) )                   
		if( do_second_scalar ) then
			s2(:,:,:)  =  s2(:,:,:)  + dt12*( 23.d0*explicit_rhs(:,:,:,5,MM0)    &
                                     - 16.d0*explicit_rhs(:,:,:,5,MM1)           &
                                     + 5.d0*explicit_rhs(:,:,:,5,MM2) )
		endif
		if(AB_order>3) step_flag='AB4'

	elseif( step_flag == 'AB4' ) then
		if(myid==0 .and. istep==3) write(0,*) '... 4th order AB step '

		ustar(:,:,:)  =  u(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,1,MM0)    &
                               - 59.d0*explicit_rhs(:,:,:,1,MM1)            &
                               + 37.d0*explicit_rhs(:,:,:,1,MM2)            &
                               - 9.d0*explicit_rhs(:,:,:,1,MM3)     )
  
		vstar(:,:,:)  =  v(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,2,MM0)    &
                               - 59.d0*explicit_rhs(:,:,:,2,MM1)            &
                               + 37.d0*explicit_rhs(:,:,:,2,MM2)            &
                               - 9.d0*explicit_rhs(:,:,:,2,MM3)     )
  
		wstar(:,:,:)  =  w(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,3,MM0)    &
                               - 59.d0*explicit_rhs(:,:,:,3,MM1)            &
                               + 37.d0*explicit_rhs(:,:,:,3,MM2)            &
                               - 9.d0*explicit_rhs(:,:,:,3,MM3)     ) 
 
		s1(:,:,:)  =  s1(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,4,MM0)  &
                                 - 59.d0*explicit_rhs(:,:,:,4,MM1)          &
                                 + 37.d0*explicit_rhs(:,:,:,4,MM2)          &
                                 - 9.d0*explicit_rhs(:,:,:,4,MM3)     )
  
		if( do_second_scalar ) then
			s2(:,:,:)  =  s2(:,:,:)  +  dt24*( 55.d0*explicit_rhs(:,:,:,5,MM0)  &
                                     - 59.d0*explicit_rhs(:,:,:,5,MM1)          &
                                     + 37.d0*explicit_rhs(:,:,:,5,MM2)          &
                                     - 9.d0*explicit_rhs(:,:,:,5,MM3)     )
		endif
	endif
 
	!-----------------------------------------------------------------
	! Handle the kappa(i) d2/dz2 ( s_bar ) terms analytically
	! since they are time independent and shouldn't be advanced
	! using an explicit method. 
	!-----------------------------------------------------------------
	! assume Laplacian dissipation and that high order derivs of ambient profiles are zero
	if( .not. high_order_operators ) then
 
 		if( ambient_profile(1) ) then
			do k=1,array_size(KDIM,YBLOCK,myid)
    			k_global = global_z_indices(START,YBLOCK,myid) + k - 1
    			s1(:,:,k) = s1(:,:,k)  + dt * kappa(1) * s1_bar(k_global,3)
   			enddo
  		endif
  	
  		if( .not. do_second_scalar ) return
  		if( ambient_profile(2) ) then
			do k=1,array_size(KDIM,YBLOCK,myid)
    			k_global = global_z_indices(START,YBLOCK,myid) + k - 1
    			s2(:,:,k) = s2(:,:,k)  + dt * kappa(2) * s2_bar(k_global,3)
   			enddo
  		endif
  	
	 endif
   
 return
end subroutine explicit_step


