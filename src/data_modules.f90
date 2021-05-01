!!=============================================
!!   Data Modules
!!=============================================

	module mpi_params
		integer,parameter               :: maxprocs=16384
		integer                         :: myid
		integer                         :: comm
		integer                         :: ierr
		integer                         :: numprocs
		integer                         :: nreqs
		integer                         :: source
		integer                         :: dest
		integer                         :: tag
	end module mpi_params

	
	module timing
		real(kind=8)                    :: t_start_time_step
		real(kind=8)                    :: t_end_time_step
		real(kind=8)                    :: t_time_step
		real(kind=8)                    :: t_total=0.d0
		real(kind=8),external           :: mpi_wtime
	end module timing

	
	module methods_params
		integer                         :: AB_order
		logical                         :: do_forcing
		logical                         :: do_second_scalar
		logical                         :: high_order_operators
		logical                         :: forcing_key(5)=.TRUE.
		logical                         :: restart=.FALSE.
		logical                         :: ambient_profile(2)=.FALSE.
		logical                         :: do_sponging=.FALSE.
	end module methods_params


	module independent_variables
		integer                         :: nx     ! global number of grid points
		integer                         :: ny     ! global number of grid points
		integer                         :: nz     ! global number of grid points
		real(kind=8)                    :: Lx     ! [m]
		real(kind=8)                    :: Ly     ! [m]
		real(kind=8)                    :: Lz     ! [m]
		real(kind=8)                    :: t_secs ! [s]
		real(kind=8)                    :: dt     ! [s]
		real(kind=8)                    :: t0     ! [s]
		real(kind=8)                    :: tf     ! [s]
		real(kind=8)                    :: tn     ! [s]
		real(kind=8)                    :: tnp1   ! [s]
		real(kind=8),allocatable        :: x(:)   ! [m] x(:), nx values
		real(kind=8),allocatable        :: y(:)   ! [m] y(:), ny values 
		real(kind=8),allocatable        :: z(:)   ! [m] z(:), nz values
		logical                         :: x_periodic
		logical                         :: y_periodic
		logical                         :: z_periodic
	contains
		subroutine initialize_coord(x,nx,Lx,periodic)
			implicit none
			integer         :: nx,i
			real(kind=8)    :: x(nx),Lx,dx
			logical         :: periodic
			if( nx > 1 ) then
				if( periodic ) then
					dx = Lx/dfloat(nx)  ! open interval[0,Lx)
				else
					dx = Lx/(nx-1.d0)   ! closed interval[0,Lx]
				endif
				do i=1,nx
					x(i) = (i-1.d0)*dx 
				enddo
			else
				x(1)=0.d0
			endif			
		end subroutine initialize_coord		
	end module independent_variables

	
	module dependent_variables
		real(kind=8),allocatable        :: u(:,:,:)        !! local portion for myid
		real(kind=8),allocatable        :: v(:,:,:)        !!      " "
		real(kind=8),allocatable        :: w(:,:,:)        !!      " "
		real(kind=8),allocatable        :: s1(:,:,:)       !!      " "
		real(kind=8),allocatable        :: s2(:,:,:)       !!      " "
		real(kind=8),allocatable        :: s1_bar(:,:)     !! global, i.e. ===>(nz,3)
		real(kind=8),allocatable        :: s2_bar(:,:)     !!      " "                       
		character(len=1)                :: scalar_kind(2)  !!  't', 's', 'p', 'r'  (passive, rho)
	end module dependent_variables
  
  
	module dimensional_scales
		real(kind=8)                    :: rho_0, rho0
		real(kind=8)                    :: nu, nu_star(3), kappa_star(3,2) ! x,y and z dirs, scalar 1/2
		real(kind=8)                    :: kappa(2)                        ! scalar 1 and 2 
		real(kind=8)                    :: T_diff(3)                       ! diffusive time scale x,y,z directions
		real(kind=8)                    :: scalar_scale(2)
		real(kind=8)                    :: pressure_scale
		real(kind=8)                    :: time_scale
		real(kind=8)                    :: velocity_scale, u0
		real(kind=8)                    :: length_scale
		real(kind=8)                    :: density_scale
		real(kind=8)                    :: f0, coriolis
		real(kind=8)                    :: g, gravity
		character(len=80)               :: scalar_name(2)		
	end module dimensional_scales


	module dimensionless_params
		integer                         :: p(3)    ! 1/2 order of diffusion operators, x,y,z
		real(kind=8)                    :: Ro
		real(kind=8)                    :: Ri
		real(kind=8)                    :: Bu
		real(kind=8)                    :: Ra
		real(kind=8)                    :: Ek
		real(kind=8)                    :: Re
		real(kind=8)                    :: Pr(2)
	end module dimensionless_params


	module differentiation_params
		integer                         :: Q=5                                         ! (Q+1)/2 term expansion in Bernoulli
		real(kind=8)                    :: filter_fraction=0.05                        ! can be 0, neg. value invokes 2/3 rule
		real(kind=8),allocatable,target :: LU_x(:,:,:), LU_y(:,:,:), LU_z(:,:,:)       ! original & factored B matrices
		integer,allocatable,target      :: ipiv_x(:,:), ipiv_y(:,:), ipiv_z(:,:)       ! corresponding pivot matrices
		real(kind=8),allocatable,target :: U_x(:,:,:), dU_x(:,:,:)                     ! U_n(x) and d/dx U_n(x)  (x,n,K)  K=1,2 for a=0,Lx
		real(kind=8),allocatable,target :: U_y(:,:,:), dU_y(:,:,:)                     ! U_n(y) and d/dy U_n(y)  (y,n,K)  K=1,2 for a=0,Ly
		real(kind=8),allocatable,target :: U_z(:,:,:), dU_z(:,:,:)                     ! U_n(z) and d/dz U_n(z)  (z,n,K)  K=1,2 for a=0,Lz
		real(kind=8),allocatable,target :: kx(:),ky(:),kz(:)                           ! wavenumber arrays, all pos. values
		real(kind=8),allocatable,target :: kxfilter(:), kyfilter(:), kzfilter(:)       ! wavenumber filters (between 0,1)
		integer(kind=8)                 :: fourier_plan(3,2), cos_plan(3), sin_plan(3) ! FFTW3 plans in x,y,z directions
		integer(kind=8)                 :: xy_plan(2,2)                                ! FFTW3 plans for xy transforms, forward and inverse
		logical                         :: fourier_done(3)=.FALSE.
		logical                         :: sin_done(3)=.FALSE.
		logical                         :: cos_done(3)=.FALSE.
	end module differentiation_params


	module etc
		character(len=80)               :: runlabel
		character(len=80)               :: logfile='output/logfile'
		character(len=80)               :: memoryfile='output/memlog'
		character(len=80)               :: message
		character(len=80)               :: step_flag='euler'
		integer                         :: istep
		integer                         :: istart
		integer                         :: iend
		integer                         :: MM0=1,MM1=2,MM2=3,MM3=4  !! for AB timestepping
		integer                         :: N=1,NM1=2                !! for AM timestepping
		logical                         :: integrate=.TRUE.
	end module etc  
  
	module intermediate_variables
		!--------------------------------------------
		!  tmp arrays in each decomposition format
		!--------------------------------------------
		real(kind=8),allocatable,dimension(:,:,:,:)    :: tmpX         ! XBLOCK format 
		real(kind=8),allocatable,dimension(:,:,:,:)    :: tmpY         ! YBLOCK format
		real(kind=8),allocatable,dimension(:,:,:,:)    :: tmpZ         ! ZBLOCK format
		real(kind=8),allocatable,dimension(:,:,:)      :: tmpYZ        ! YZ plane of YBLOCK format
		!--------------------------------------------
		! 3d arrays in YBLOCK format
		!--------------------------------------------
		real(kind=8),allocatable,dimension(:,:,:)      :: phi          ! YBLOCK format
		real(kind=8),allocatable,dimension(:,:,:)      :: ustar        ! YBLOCK format
		real(kind=8),allocatable,dimension(:,:,:)      :: vstar        ! YBLOCK format
		real(kind=8),allocatable,dimension(:,:,:)      :: wstar        ! YBLOCK format
		real(kind=8),allocatable,dimension(:,:,:)      :: div_u        ! YBLOCK format
		real(kind=8),allocatable,dimension(:,:,:,:,:)  :: explicit_rhs ! YBLOCK format
		real(kind=8),allocatable,dimension(:,:,:,:,:)  :: implicit_rhs ! YBLOCK format 
	end module intermediate_variables



	module decomposition_params
		integer                   :: np   ! np=p1*p2
		integer                   :: p1
		integer                   :: p2
		!------------------------------------------------------------------------------------------
		! mem_order(1,XBLOCK) = the fortran array index corresponding to x in XBLOCK decomposition
		! mem_order(2,XBLOCK) = the fortran array index corresponding to y in XBLOCK decomposition
		! mem_order(3,XBLOCK) = the fortran array index corresponding to z in XBLOCK decomposition
		!------------------------------------------------------------------------------------------
		integer                   :: mem_order(3,3)
		!-----------------------------------------------------------------------------------------
		!  e.g. layout(XCOORD,ZBLOCK) defines the number of processors used in the decomposition of
		!                             x when the global data set is arranged in ZBLOCKS
		!-----------------------------------------------------------------------------------------
		integer                   :: layout(3,3)
		!-----------------------------------------------------------------------------------------
		!  e.g. proc_row(YBLOCK,pid) defines the row assignment of pid in the n1xn2 processor grid
		!                            when the global data set is arranged in YBLOCKS
		!-----------------------------------------------------------------------------------------
		integer,allocatable       :: proc_row(:,:)
		!-----------------------------------------------------------------------------------------
		!  e.g. proc_col(YBLOCK,pid) defines the col assignment of pid in the n1xn2 processor grid
		!                            when the global data set is arranged in YBLOCKS
		!-----------------------------------------------------------------------------------------
		integer,allocatable       :: proc_col(:,:)
		!-----------------------------------------------------------------------------------------
		!  e.g. array_size(JDIM,XBLOCK,pid) defines the number of grid points in the local data array
		!                                   on pid in the JDIM storage index in the XBLOCK decomposition
		!-----------------------------------------------------------------------------------------
		integer,allocatable       :: array_size(:,:,:)
		!-----------------------------------------------------------------------------------------
		!  e.g. global_x_indices(END,XBLOCK,pid) defines the the LAST global x index for the local
		!                                   data array on pid in the XBLOCK decomposition
		!-----------------------------------------------------------------------------------------
		integer,allocatable       :: global_x_indices(:,:,:)
		integer,allocatable       :: global_y_indices(:,:,:)
		integer,allocatable       :: global_z_indices(:,:,:)
		!-----------------------------------------------------------------------------------------
		!  these are simply conventions and should never change
		!  (to keep me from having to use 1,2,3 in different contexts)
		!-----------------------------------------------------------------------------------------
		integer,parameter         :: XCOORD=1
		integer,parameter         :: YCOORD=2
		integer,parameter         :: ZCOORD=3
		integer,parameter         :: XBLOCK=1
		integer,parameter         :: YBLOCK=2
		integer,parameter         :: ZBLOCK=3
		integer,parameter         :: IDIM=1   ! convention: fortran arrays are indexed (IDIM,JDIM,KDIM)
		integer,parameter         :: JDIM=2   ! i.e. IDIM is fastest varying array index, then JDIM,KDIM
		integer,parameter         :: KDIM=3
		integer,parameter         :: START=1
		integer,parameter         :: END=2

		integer,allocatable       :: jk_indices_YBLOCK(:,:)
		integer,allocatable       :: ijk_indices_YBLOCK(:,:)
		integer,allocatable       :: ijk_indices_ZBLOCK(:,:)
  contains
  		subroutine get_my_xvals( xvals, iblock, myid )
  			use independent_variables,     only: x
  			implicit none
  			integer, intent(in)               :: iblock, myid
  			integer                           :: nvals, i, ig
  			real(kind=8), intent(out)         :: xvals(:)  			                             
			nvals = global_x_indices(END  ,iblock,myid) \
			      - global_x_indices(START,iblock,myid) + 1			
			do i = 1,nvals
				ig = global_x_indices(START,iblock,myid) + i - 1
				xvals(i) = x(ig)
			enddo			
		end subroutine get_my_xvals
		
		subroutine get_my_yvals( yvals, iblock, myid )
  			use independent_variables,     only: y
  			implicit none
  			integer, intent(in)               :: iblock, myid
  			integer                           :: nvals, i, ig
  			real(kind=8), intent(out)         :: yvals(:)  			                             
			nvals = global_y_indices(END  ,iblock,myid) \
			      - global_y_indices(START,iblock,myid) + 1			
			do i = 1,nvals
				ig = global_y_indices(START,iblock,myid) + i - 1
				yvals(i) = y(ig)
			enddo			
		end subroutine get_my_yvals
		
		subroutine get_my_zvals( zvals, iblock, myid )
  			use independent_variables,     only: z
  			implicit none
  			integer, intent(in)               :: iblock, myid
  			integer                           :: nvals, i, ig
  			real(kind=8), intent(out)         :: zvals(:)  			                             
			nvals = global_z_indices(END  ,iblock,myid) \
			      - global_z_indices(START,iblock,myid) + 1			
			do i = 1,nvals
				ig = global_z_indices(START,iblock,myid) + i - 1
				zvals(i) = z(ig)
			enddo			
		end subroutine get_my_zvals
  end module decomposition_params




	module io_params
		! data read in from input/io_params
		integer,parameter             :: maxsets=2048
		integer                       :: num_file_sets             ! user supplied
		character(len=80)             :: filename_root(maxsets)    ! user supplied
		character(len=80)             :: mode(maxsets)             ! user supplied
		integer                       :: ilocs(3,maxsets)          ! user supplied
		integer                       :: jlocs(3,maxsets)          ! user supplied
		integer                       :: klocs(3,maxsets)          ! user supplied
		integer                       :: variable_key(11,maxsets)  ! user supplied
		integer                       :: nsteps(maxsets)           ! user supplied

		! nonsingleton indices arranged 1st in the contiguous array 'outdata' that
		! contains only and all the data to be actually written to the netcdf file
		integer                       :: dimid_x(maxsets)  ! 'outdata' array dimension for x in nc file
		integer                       :: dimid_y(maxsets)  ! 'outdata' array dimension for x in nc file
		integer                       :: dimid_z(maxsets)  ! 'outdata' array dimension for x in nc file
		integer                       :: dimid_t(maxsets)  ! 'outdata' array dimension for x in nc file

		! these arrays are used in the actual netcdf write statements
		! to define how to interpret the layout of the 'outdata' array
		integer                       :: nspace(maxsets)        ! number of nonsingleton space dimensions
		integer                       :: time_counter(maxsets)  ! time index
		integer                       :: count(4,maxsets)

		! filename constructed given its root, num of dimensions and topdir
		character(len=80)             :: fullname(maxsets)

		! indices, counts of local YBLOCK data to extract for writing
		integer                       :: my_x0(maxsets)    ! initial local index for x coordinate
		integer                       :: my_x1(maxsets)    ! final   local index for x coordinate
		integer                       :: my_xinc(maxsets)  ! increment for x coordinate sampling
		integer                       :: my_nx(maxsets)    ! number of x coord vals to be written

		integer                       :: my_y0(maxsets)    ! initial local indey for y coordinate
		integer                       :: my_y1(maxsets)    ! final   local indey for y coordinate
		integer                       :: my_yinc(maxsets)  ! increment for y coordinate sampling
		integer                       :: my_ny(maxsets)    ! number of y coord vals to be written

		integer                       :: my_z0(maxsets)    ! initial local indez for z coordinate
		integer                       :: my_z1(maxsets)    ! final   local indez for z coordinate
		integer                       :: my_zinc(maxsets)  ! increment for z coordinate sampling
		integer                       :: my_nz(maxsets)    ! number of z coord vals to be written

		! key specifying whether there is data to be written
		logical                       :: do_write(maxsets)
		logical                       :: write_s1_bar(maxsets)
		logical                       :: write_s2_bar(maxsets)
	end module io_params


	module boundary_data
		real(kind=8), allocatable    ::   east_vals(:,:,:)      ! (y,z,id) at x=0
		real(kind=8), allocatable    ::   west_vals(:,:,:)      ! (y,z,id) at x=Lx
		real(kind=8), allocatable    ::  south_vals(:,:,:)      ! (x,z,id) at y=0
		real(kind=8), allocatable    ::  north_vals(:,:,:)      ! (x,z,id) at y=Ly
		real(kind=8), allocatable    :: bottom_vals(:,:,:)      ! (x,y,id) at z=0
		real(kind=8), allocatable    ::    top_vals(:,:,:)      ! (x,y,id) at z=Lz
		
		real(kind=8), allocatable    ::   east_derivs(:,:,:)    ! (y,z,id) at x=0
		real(kind=8), allocatable    ::   west_derivs(:,:,:)    ! (y,z,id) at x=Lx
		real(kind=8), allocatable    ::  south_derivs(:,:,:)    ! (x,z,id) at y=0
		real(kind=8), allocatable    ::  north_derivs(:,:,:)    ! (x,z,id) at y=Ly
		real(kind=8), allocatable    :: bottom_derivs(:,:,:)    ! (x,y,id) at z=0
		real(kind=8), allocatable    ::    top_derivs(:,:,:)    ! (x,y,id) at z=Lz	
		real(kind=8), allocatable    :: step_e(:),step_w(:)     ! (x)
		real(kind=8), allocatable    :: step_s(:),step_n(:)     ! (y)
		real(kind=8), allocatable    :: step_b(:),step_t(:)     ! (y)
	end module boundary_data

