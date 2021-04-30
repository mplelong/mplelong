subroutine InitialConditions  
	use etc,                    only: logfile
	use independent_variables,  only: x,y,z
	use methods_params,         only: do_second_scalar,ambient_profile
	use decomposition_params 
	use intermediate_variables
	use dependent_variables
	use dimensional_scales
	use mpi_params
	implicit none 
 
	real(kind=8),allocatable   :: my_xvals(:)
	real(kind=8),allocatable   :: my_yvals(:)
	real(kind=8),allocatable   :: my_zvals(:)
	integer                    :: id,m1,m2,m3,i,j,k 
 
  
	if(myid==0) then
		write(0,*) ' ................'
		write(0,*) ' ................     hello world from InitialConditions'
		open(1,file=logfile,position='append') 
		write(1,*) '  '
		write(1,*) '  '
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '                  InitialConditions Report:'
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '  '
		close(1)
	endif
  		 
	!------------------------------------------------
	!    YBLOCK arranged as data(y,x,z)
	!    data(m2,m1,m3)
	!------------------------------------------------
	m1=array_size(JDIM,YBLOCK,myid)   ! num of x data points in YBLOCK format (index=2=JDIM)
	m2=array_size(IDIM,YBLOCK,myid)   ! num of y data points in YBLOCK format (index=1=IDIM)
	m3=array_size(KDIM,YBLOCK,myid)   ! num of z data points in YBLOCK format
	allocate( my_xvals(m1) ) ; my_xvals(:)=0.d0
	allocate( my_yvals(m2) ) ; my_yvals(:)=0.d0
	allocate( my_zvals(m3) ) ; my_zvals(:)=0.d0
	
	call get_my_xvals( my_xvals, YBLOCK, myid )
	call get_my_yvals( my_yvals, YBLOCK, myid )
	call get_my_zvals( my_zvals, YBLOCK, myid )
 
  		 
 	id=1
 	call user_ics(my_xvals,my_yvals,my_zvals,u,id,m1,m2,m3)
 
 
 	id=2
	call user_ics(my_xvals,my_yvals,my_zvals,v,id,m1,m2,m3)

 
	id=3
 	call user_ics(my_xvals,my_yvals,my_zvals,w,id,m1,m2,m3)
 
  
	id=4
	call user_ics(my_xvals,my_yvals,my_zvals,s1,id,m1,m2,m3) 
 
 
	if( do_second_scalar ) then 
  		id=5
  		call user_ics(my_xvals,my_yvals,my_zvals,s2,id,m1,m2,m3)
 	endif
 
	deallocate( my_xvals )
	deallocate( my_yvals )
	deallocate( my_zvals ) 
 
 
 	if(myid==0) then
  		open(1,file=logfile,position='append')
  			if( do_second_scalar ) then
  				write(1,*) ' ................      initialization of YBLOCK arrays u,v,w,s1,s2:  '
  			else
  				write(1,*) ' ................      initialization of YBLOCK arrays u,v,w,s1:     '
  			endif
  			write(1,*) ' ................      initialization via flow_solve_user_routines.f90/user_ics '
  			write(1,*) ' -----> InitialConditions routine exiting normally  <---------- '
  		close(1)
	endif
 
	call mpi_barrier(comm,ierr)
 return
end subroutine InitialConditions
