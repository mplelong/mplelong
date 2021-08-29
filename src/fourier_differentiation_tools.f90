!
!----------------------------------------------------------------------
! Module to use FFTW3 to differentiate 1d data through term by term
! differentiation of the appropriate trigonometric series expansion.
!
! Written by Kraig Winters, Scripps Institution of Oceanography, UCSD
!----------------------------------------------------------------------
!
!-------------------------------------------------------------------------
   MODULE fourier_differentiation_tools
!-------------------------------------------------------------------------
	use differentiation_params
	use independent_variables,  only: x,y,z,nx,ny,nz,Lx,Ly,Lz  ! (global coord. arrays)
    implicit none

  CONTAINS
  
!-------------------------------------------------------------------------
SUBROUTINE get_kind(in,out,type,fft_kind,label)
!-------------------------------------------------------------------------
!
!   determine the "kind" of FFTW3  transform we'll need
!
!-------------------------------------------------------------------------
	implicit none
	include 'fftw3.f'

	integer,intent(out)               :: fft_kind
 	character (len=80),intent(out)    :: label
 	character (len=11),intent(in)     :: in
 	character (len=11),intent(in)     :: out
 	character (len=11),intent(in)     :: type
 	
 	if( trim(type)=='fourier' ) then
		if( in=='real' .and. out=='complex') then
			fft_kind=FFTW_R2HC
			label='no flag needed'
		elseif( in=='complex' .and. out=='real') then
			fft_kind=FFTW_HC2R
			label='no flag needed'
		elseif( in=='real' .and. out=='halfcomplex') then
			fft_kind=FFTW_R2HC
			label='FFTW_R2HC'
		elseif( in=='halfcomplex' .and. out=='real') then
			fft_kind=FFTW_HC2R
			label='FFTW_HC2R'
		elseif( in=='complex' .and. out=='complex') then
			fft_kind=-999
			label='no flag needed'
  		else
			write(0,*) 'problem getting kind of fourier transform',in,out
			stop
		endif        
	endif
	
       
 	if( trim(type)=='sin' ) then
		if( in=='real' .and. out=='real') then
			fft_kind=FFTW_RODFT00
			label='FFTW_RODFT00'
		elseif( in=='complex' .and. out=='complex') then
			fft_kind=-999
			label='no label needed'
		else
			write(0,*) 'problem getting kind of fourier transform',in,out
			stop
		endif  
	endif
	
 
	if( trim(type)=='cos' ) then
		if( in=='real' .and. out=='real') then
			fft_kind=FFTW_REDFT00
			label='FFTW_REDFT00'
		elseif( in=='complex' .and. out=='complex') then
			fft_kind=-999
			label='no label needed'
		else
			write(0,*) 'problem getting kind of fourier transform',in,out
			stop
		endif  
	endif
  	   	
return
END SUBROUTINE get_kind



!-------------------------------------------------------------------------
SUBROUTINE fourier_init(exp_type,n,Lx,kx,kfilter,plan_f,plan_i)
!-------------------------------------------------------------------------
!
!  initialize wavenumber & filter arrays, forward and inverse FFTW3 plans
!
!-------------------------------------------------------------------------
	use mpi_params,                  only: myid
	use differentiation_params,      only: filter_fraction
	implicit none
	
	real(kind=8), external              :: myexp    ! decaying exponential w/o underflow
	character(len=80),intent(in)        :: exp_type
	integer, intent(in)                 :: n
	real(kind=8), intent(in)            :: Lx
	real(kind=8), intent(out)           :: kx(n),kfilter(n)
	integer(kind=8), intent(out)        :: plan_f,plan_i
	
	real(kind=8)                        :: pi,dk,kmax,gamma
	real(kind=8),allocatable            :: in(:), out(:)
	integer                             :: fft_kind
	character(len=80)                   :: reality_in, reality_out, label, this_type
	integer                             :: i,k
	integer                             :: rank_trans
	integer                             :: howmany_rank
	integer                             :: n_trans,   n_loop
	integer                             :: is_trans, is_loop
	integer                             :: os_trans, os_loop
	integer(kind=8)                     :: plan, plan_sin, plan_cos
	include 'fftw3.f'
	
	if( n < 4 ) then
		if(myid==0) then
			write(0,*) ' ................      Warning: fourier initialization for n < 4: ',  &
			trim( exp_type)
		endif
		kx=0.d0
		kfilter=1.d0
		plan_f=-999
		plan_i=-999
		return  
	endif
     
     
	if( trim(exp_type)=='fourier' .and. mod(n,2) /= 0 ) then
		write(0,*) 'n = ',n
		stop 'need even number of gridpoints w/ 1 end excluded for fourier implementation'
	elseif( trim(exp_type)=='cos' .and. mod(n,2) /= 1 ) then
		stop 'need odd number of gridpoints w/ both ends included for cos implementation'
	elseif( trim(exp_type)=='sin' .and. mod(n,2) /= 1 ) then
		stop 'need odd number of gridpoints w/ both ends included for sin implementation'
	endif
    
	rank_trans=1      !! 1d transforms
	howmany_rank=1    !! Any loops over additional dimensions 
                      !! specified as though data laid out
                      !! in 1d array in memory (not an issue here)
  
	pi=4.d0*atan(1.d0)
	
	if( trim(exp_type) .eq. 'sin' .or. trim(exp_type) .eq. 'cos') then
  
  		!  initialize wavenumber and wavenumber filter arrays
    	dk = pi/Lx
    	kmax = (n-1.)*dk
    	gamma = filter_fraction*n*dk
    	do i=1,n
    		kx(i) = (i-1.)*dk
    		if( filter_fraction > 0.d0 ) then
    			kfilter(i) = 1.d0 - myexp(-((kx(i)-kmax)/gamma)**2)
    		elseif( filter_fraction < 0.d0 ) then
    			if( kx(i) > (2.d0/3.d0)*kmax ) then
    				kfilter(i) = 0.d0
    			else
    				kfilter(i) = 1.d0
    			endif
    		elseif( filter_fraction == 0.d0 ) then
    			kfilter(i) = 1.d0
    		endif
    	enddo
    
		!!=======================================================
		!!make a plan for a sin transform ignoring zero endvalues
		!!=======================================================
    
		!!describe the transforms  in dim 1
		n_trans=n-2              ! length of transform excluding zeros at ends
		is_trans=1               ! stride btwn elements in a single transform
		os_trans=is_trans        ! set output stride = input stride
        
		!!describe the loops over additional dims, i.e. 2,3
		n_loop=1                ! collapse more general scheme to 1d 
		is_loop=1               ! stride betwn transforms
		os_loop=is_loop         ! set output stride = input stride
    
		reality_in=  'real'
		reality_out= 'real'
		this_type='sin'
    	call get_kind(reality_in,reality_out,this_type,fft_kind,label)
    
		allocate( in(n),out(n) )
		in = 0.d0
		out = 0.d0
		call dfftw_plan_guru_r2r(plan,             &
                             	rank_trans,        &
                             	n_trans,           &
                             	is_trans,          &
                             	os_trans,          &
                             	howmany_rank,      &
                             	n_loop,            &
                             	is_loop,           &
                             	os_loop,           &
                             	in(2),             &
                             	out(2),            &
                             	fft_kind,          &
                             	FFTW_EXHAUSTIVE)
		if( trim(exp_type) == 'sin' ) plan_f=plan
		if( trim(exp_type) == 'cos' ) plan_i=plan    
    	plan_sin = plan_f   ! sin w/ no endpoints, testing only
    
		if(myid==0) then
			!write(0,*) ' ................       testing ',plan,trim(exp_type),n
			!call dfftw_execute_r2r(plan,in(2),out(2))
		endif 
		deallocate(in,out)
    
    
		!!=======================================================
		!! Now make a plan for a cos transform keeping endvalues
		!!=======================================================
    
		!!describe the transforms  in dim 1
		n_trans=n                ! length of transform including zeros at ends
		is_trans=1               ! stride btwn elements in a single transform
		os_trans=is_trans        ! set output stride = input stride
        
		!!describe the loops over additional dims, i.e. 2,3
		n_loop=1                ! collapse more general scheme to 1d 
		is_loop=1               ! stride betwn transforms
		os_loop=is_loop         ! set output stride = input stride
    
		reality_in=  'real'
		reality_out= 'real'
		this_type='cos'
		call get_kind(reality_in,reality_out,this_type,fft_kind,label)

		allocate( in(n),out(n) )
		in = 0.d0
		out =0.d0
		call dfftw_plan_guru_r2r(plan,              &
                             	rank_trans,         &
                             	n_trans,            &
                             	is_trans,           &
                            	os_trans,           &
                            	howmany_rank,       &
                            	n_loop,             &
                            	is_loop,            &
                            	os_loop,            &
                            	in,                 &
                             	out,                &
                             	fft_kind,           &
                             	FFTW_EXHAUSTIVE)
		if( trim(exp_type) == 'cos' ) plan_f=plan
		if( trim(exp_type) == 'sin' ) plan_i=plan    
    	plan_cos = plan_f   ! testing only
    
		! check transform/inverse transform pair for sin
		! do i=1,n
		!  in(i) = sin(2.*pi*(i-1.)/dfloat(n-1.))
		! enddo
    
		!write(0,*) plan_sin,plan_cos
		if(myid==0) then
    		!write(0,*) ' ................       testing ',plan,trim(exp_type),n
			!call dfftw_execute_r2r(plan_sin,in(2),out(2))
			!call dfftw_execute_r2r(plan_cos,out(1),in(1))
		endif    
		!do i=1,n
			! write(0,*) i,sin(2.*pi*(i-1.)/dfloat(n-1.)),in(i)/(2.*n)
		!enddo   
    	deallocate(in,out)
     
	elseif( trim(exp_type) == 'fourier' ) then
  
		dk = 2.d0*pi/Lx   
    	kmax = (n/2+1.d0)*dk
    	gamma = filter_fraction*n*dk
		do i=1,n/2+1       
			kx(i) = (i-1.)*dk
    		if( filter_fraction > 0.d0 ) then
    			kfilter(i) = 1.d0 - myexp(-((kx(i)-kmax)/gamma)**2)
    		elseif( filter_fraction < 0.d0 ) then
    			if( kx(i) > (2.d0/3.d0)*kmax ) then
    				kfilter(i) = 0.d0
    			else
    				kfilter(i) = 1.d0
    			endif
    		elseif( filter_fraction == 0.d0 ) then
    			kfilter(i) = 1.d0
    		endif
		enddo
     
		do i=n,n/2+2,-1
			kx(i) = -kx( n-i+2 )     
			kfilter(i) = kfilter(n-i+2)
		enddo
   
		!!describe the transforms  in dim 1
		n_trans=n                ! length of transform
		is_trans=1               ! stride btwn elements in a single transform
		os_trans=is_trans        ! set output stride = input stride
        
		!!describe the loops over additional dims, i.e. 2,3
		n_loop=1                ! collapse more general scheme to 1d 
		is_loop=1               ! stride betwn transforms
		os_loop=is_loop         ! set output stride = input stride
    
		!! construct the forward plan
		reality_in='real'
		reality_out='halfcomplex'
		call get_kind(reality_in,reality_out,exp_type,fft_kind,label)
		allocate( in(n),out(n) )
		in = 0.d0
		out =0.d0
		call dfftw_plan_guru_r2r(plan_f,        &
                             rank_trans,        &
                             n_trans,           &
                             is_trans,          &
                             os_trans,          &
                             howmany_rank,      &
                             n_loop,            &
                             is_loop,           &
                             os_loop,           &
                             in,                &
                             out,               &
                             fft_kind,          &
                             FFTW_EXHAUSTIVE)
		if(myid==0) then
			!write(0,*) ' ................       testing ',plan_f,trim(exp_type),n
			!call dfftw_execute_r2r(plan_f,in,out)
		endif  
    
		!! construct the inverse plan
		reality_in='halfcomplex'
		reality_out='real'
		call get_kind(reality_in,reality_out,exp_type,fft_kind,label)     
		call dfftw_plan_guru_r2r(plan_i,        &
                             rank_trans,        &
                             n_trans,           &
                             is_trans,          &
                             os_trans,          &
                             howmany_rank,      &
                             n_loop,            &
                             is_loop,           &
                             os_loop,           &
                             in,                &
                             out,               &
                             fft_kind,          &
                             FFTW_EXHAUSTIVE)
		if(myid==0) then
			!write(0,*) ' ................       testing ',plan_i,trim(exp_type),n
			!call dfftw_execute_r2r(plan_i,in,out)
		endif
		deallocate(in,out)
	endif

return
END SUBROUTINE fourier_init


!-------------------------------------------------------------------------
SUBROUTINE fourier_deriv(f,df,n,order,exp_type,kx,kfilter,tmp,plans)
!-------------------------------------------------------------------------
!
!  differentiate f(:) "order" times and store the result in df(:)
!  using fourier/cos/sin expansion & term by term differentiation 
!  
!   n      number of grid points
!   order  which order deriv to take  (not the order of convergence)
!   kx     input:   wavenumber vector
!   kfilt  input:  wavenumber filter
!   tmp    input:   workspace real(kind=8) dimension(n)
!   exp_type  'fourier','cos' or 'sin' , how f is expanded
!-------------------------------------------------------------------------
	use mpi_params,                      only: myid,comm,ierr
	implicit none
 	integer, intent(in)                     :: n,order
 	real(kind=8)                            :: f(n)
 	real(kind=8)                            :: df(n)
 	integer(kind=8)                         :: plans(2)
 	character(len=80)                       :: exp_type
 	real(kind=8)                            :: kx(n),kfilter(n)
 	real(kind=8)                            :: tmp(n)
 	
 	integer                                 :: i,istart
 	integer(kind=8)                         :: plan_f, plan_i
 	real(kind=8)                            :: normfactor,xx
	integer,parameter                       :: maxorder=12
	integer,save                            :: sign_fourier(maxorder)
	integer,save                            :: sign_cos(maxorder)
	integer,save                            :: sign_sin(maxorder)
	logical,save                            :: do_phase_shift(maxorder)
	logical,save                            :: first_entry=.TRUE.

	if(n<4) stop 'Call to fourier_deriv with vector of length < 4'
	
	if( first_entry ) then
		!! sign and phase shift params, fourier method  
		sign_fourier(1:4) = (/-1,-1,1,1/)
		sign_fourier(5:8) = (/-1,-1,1,1/)
		sign_fourier(9:12) = (/-1,-1,1,1/)
		do_phase_shift(1:2:maxorder-1)=.TRUE.
		do_phase_shift(2:2:maxorder ) =.FALSE.
		!! sign parameters cos method, input data is assumed even
		sign_cos(1:4) = (/-1,-1,1,1/)     ! cos:  --> -sin --> -cos --> sin --> cos
		sign_cos(5:8) = (/-1,-1,1,1/)
		sign_cos(9:12) = (/-1,-1,1,1/) 
		!! sign parameters sin method, input data is assumed odd
		sign_sin(1:4) = (/1,-1,-1,1/)     ! sin:  --> cos --> -sin --> -cos --> sin
		sign_sin(5:8) = (/1,-1,-1,1/)
		sign_sin(9:12) = (/1,-1,-1,1/)
		first_entry=.FALSE.
	endif
		

	if( trim(exp_type) == 'fourier' ) then
 
		normfactor = (1.d0/dfloat(n)) * sign_fourier(order)
		plan_f = plans(1)
		plan_i = plans(2)
  
		!-----------------------------------------------
		! do the forward transform 
		!-----------------------------------------------
		call dfftw_execute_r2r(plan_f,f,tmp)  !! fourier
 
		!-------------------------------------------------------------
		! if necessary, the "i" part of -ikx
		! location shift necessary for real--> 1/2 complex transform
		! implementation, other ways also possible
		!-------------------------------------------------------------
		if( do_phase_shift(order) ) then  
			do i=2,n/2
				xx = tmp(i)
				tmp(i) = tmp(n-i+2)
				tmp(n-i+2) = xx
			enddo
		endif
  
		!-------------------------------------------------------------
		! now the "-k_x" and the normalization:
		! tmp(:) = - kx(:)*normfactor*tmp(:)
		! ---> generalized for arbitrary orders, see above
		! ---> wavenumber filtering applied
		!-------------------------------------------------------------
		do i=1,n
			tmp(i) = kfilter(i)*tmp(i)*normfactor*(kx(i))**order
		enddo
		tmp(n/2+1) = 0.d0    !! zero the nyquist frequency
 
		!-------------------------------------------------------------
		! do the inverse transform
		!-------------------------------------------------------------
		call dfftw_execute_r2r(plan_i,tmp,df)

	elseif( trim(exp_type) == 'cos' ) then
    
		normfactor = (1.d0/(2.d0*(dfloat(n)-1.d0))) * sign_cos(order)
		plan_f = plans(1)
		if( mod(order,2)==0 ) then   ! inverse is also n point cosine transform
			plan_i = plan_f
			istart = 1
		else                         ! inverse is shortened sin transform
			plan_i = plans(2)
			istart = 2             
		endif
  
		call dfftw_execute_r2r(plan_f,f,tmp)        !! forward/cos transform
  
		!-----------------------------------------------
		! wavenumber multiplication and filtering
		!-----------------------------------------------
		do i=1,n
			tmp(i) = kfilter(i)*tmp(i)*normfactor*(kx(i))**order
		enddo
		tmp(n)=0.d0      ! set nyquist to zero
  
		!-------------------------------------------------------------
		! do the inverse transform
		!-------------------------------------------------------------
		call dfftw_execute_r2r(plan_i,tmp(istart),df(istart))
  
		!-------------------------------------------------------------
		! if final result is sin expandable, i.e. odd deriv of even fn
		!  ==> zero the endvals
		!-------------------------------------------------------------
		if( order==1 .or. mod(order,2) .ne. 0 ) then  
			df(1)=0.d0
			df(n)=0.d0
		endif
  
	elseif( trim(exp_type) == 'sin' ) then
    
		normfactor = 1.d0/(2.d0*(dfloat(n)-1.d0)) * sign_sin(order)
		plan_f=plans(1)
		if( mod(order,2)==0 ) then  ! inverse is also shortened sin transform
			plan_i = plan_f
			istart = 2
		else                        ! inverse is n point cosine transform
			plan_i = plans(2)
			istart = 1
		endif
  
		!-----------------------------------------------
		!  forward/sin transform
		!-----------------------------------------------
		call dfftw_execute_r2r(plan_f,f(2),tmp(2))
  
		!-----------------------------------------------
		! wavenumber multiplication and filtering
		!-----------------------------------------------
		do i=1,n
			tmp(i) = kfilter(i)*tmp(i)*normfactor*(kx(i))**order
		enddo
		tmp(n)=0.d0    !! nyquist
  
		!-------------------------------------------------------------
		! do the inverse transform
		!-------------------------------------------------------------
		call dfftw_execute_r2r(plan_i,tmp(istart),df(istart))
  
		!-------------------------------------------------------------
		! if final result is sin expandable, i.e. even deriv of odd fn
		!  ==> zero the endvals
		!-------------------------------------------------------------
		if( mod(order,2) == 0) then
			df(1) = 0.d0
			df(n) = 0.d0
		endif
  
	endif
	  
return
END SUBROUTINE fourier_deriv


!-------------------------------------------------------------------------
SUBROUTINE differentiate_fcs(f,df,n,dir,method,order)
!-------------------------------------------------------------------------
! compute the "order" derivative of a 1d array of data using term by term 
! differentiation of the F/C/S expansion
! method is the expansion type for the data 'fourier', 'cos' or 'sin'
!-------------------------------------------------------------------------
	use independent_variables,        only: nx,ny,nz
	use differentiation_params,       only: kx,ky,kz,kxfilter,kyfilter,kzfilter
	use differentiation_params,       only: cos_plan,sin_plan,fourier_plan
	implicit none 
	integer, intent(in)                  :: n
	integer, intent(in)                  :: dir       ! 1,2,3 for diff wrt to x,y,z coordinate
	integer, intent(in)                  :: order
	real(kind=8)                         :: f(n)   
	real(kind=8)                         :: df(n)
	character(len=80), intent(in)        :: method   ! method of 1st deriv, even when order>1 	 
	integer(kind=8)                      :: plans(2)
	character(len=80)                    :: exp_type
	integer                              :: nmax
	real(kind=8),allocatable,save        :: k(:),kfilter(:),tmp(:)
	logical, save                        :: first_entry=.TRUE.
	
	if(first_entry) then		
		nmax = MAXVAL( (/nx,ny,nz/) )
		allocate( k(nmax), kfilter(nmax), tmp(nmax) )
		k = 0.d0 ; kfilter=0.d0 ; tmp = 0.d0		
		first_entry=.FALSE.
	endif
	
	if( n .ne. nx .and. n .ne. ny .and. n .ne. nz )    stop ' bad n value passed to differentiate_fcs '
	if( dir .ne. 1 .and. dir .ne. 2 .and. dir .ne. 3 ) stop ' bad dir value passed to differentiate_fcs '
	
   
	if(dir==1) then
		k(1:n)       = kx(1:nx)         ! these lengths should match
		kfilter(1:n) = kxfilter(1:nx)
	elseif(dir==2) then
		k(1:n)       = ky(1:ny)         ! these lengths should match
		kfilter(1:n) = kyfilter(1:ny)
	elseif(dir==3) then
		k(1:n)       = kz(1:nz)         ! these lengths should match
		kfilter(1:n) = kzfilter(1:nz)
	endif

	if( n==1 ) then
		df(1)=0.d0
		return
	endif
  
	exp_type = method
	if( method=='cos') then
		plans(1) = cos_plan(dir)
		plans(2) = sin_plan(dir)
	elseif( method=='sin') then
		plans(1) = sin_plan(dir)
		plans(2) = cos_plan(dir)
	elseif( method=='fourier' ) then
		plans(1) = fourier_plan(dir,1)
		plans(2) = fourier_plan(dir,2)
	endif
  	
	call fourier_deriv(f,df,n,order,exp_type,k,kfilter,tmp,plans)
	
	   
end subroutine differentiate_fcs

!-------------------------------------------------------------------------
   END MODULE fourier_differentiation_tools
!-------------------------------------------------------------------------



