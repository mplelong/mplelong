subroutine user_analysis
	use mpi_params,             only: myid
	use decomposition_params
	use dependent_variables,    only: u,v,w
	use methods_params,         only: do_second_scalar
	use independent_variables,  only: x,y,z,t_secs,tf,t0,dt
	use etc,                    only: istep,istart,iend
	use user_params
  
	implicit none
	integer,save                    ::  locnx,locny,locnz
	logical                         ::  first_entry=.TRUE. 
	include 'mpif.h'
   
	if( first_entry ) then
		locnx = array_size(JDIM,YBLOCK,myid)
		locny = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		first_entry = .FALSE.
	endif
	  
 return  
end subroutine user_analysis
