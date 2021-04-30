subroutine user_analysis
	use mpi_params,             only: myid,comm,ierr,numprocs
	use decomposition_params,   only: proc_row,proc_col,YBLOCK,array_size,IDIM,JDIM,KDIM,p1,p2
	use dependent_variables,    only: u,v,w,s1,s2
	use methods_params,         only: do_second_scalar
	use user_params,            only: X_means,rho_bar
  
	implicit none
	integer                         ::  i,j,k,id,pid
	integer,save                    ::  locnx,locny,locnz,nvars,count
	integer,save                    ::  group_world,rows_group(0:4096),members(4096)
	integer,save                    ::  rows_comm(0:4096)
	real(kind=8),allocatable        ::  tmp(:)
	logical                         ::  first_entry=.TRUE. 
	include 'mpif.h'
 
	if( array_size(JDIM,YBLOCK,myid) == 1 ) return  ! no need to x avg if nx=1
  
	if( first_entry ) then
		locnx = array_size(JDIM,YBLOCK,myid)
		locny = array_size(IDIM,YBLOCK,myid)
		locnz = array_size(KDIM,YBLOCK,myid)
		count = locny*locnz
		allocate( X_means(locny,locnz,6) )       ! f(y,z) for u,v,w,s1,s2 + scratch
		allocate( rho_bar(locnz),tmp(locnz) )    
		rho_bar(:) = 0.d0
		tmp(:) = 0.d0
   
   		nvars = 4
   		if( do_second_scalar ) nvars = 5
   
   		!---------------------------------------------
   		!  create p2 new groups for communicating
   		!  between proc_rows with the same proc_col
   		!---------------------------------------------
   		call MPI_Comm_group(comm, group_world, ierr)   ! handle to group of ALL processors
   
   		do i=0,p2-1                   ! loop through all proc_col indices
   
    		j=1                          ! get the members in new subgroup
    		do pid = 0,numprocs-1
     			if( proc_col(YBLOCK,pid) == i ) then
      				members(j) = pid           ! processor pid has proc_col equal to i 
      				j = j + 1
     			endif
    		enddo
    
    		!--------------------------------------------------------
    		! create the mpi group with proc_col equal to i
    		!--------------------------------------------------------
    		call MPI_GROUP_INCL(group_world, p1, members(1:p1), rows_group(i), ierr)
    
    		!--------------------------------------------------------
    		! create the mpi communicator with proc_col equal to i
    		!--------------------------------------------------------
    		call MPI_COMM_CREATE(comm, rows_group(i), rows_comm(i),  ierr)
    
   		enddo
   
   		!----------------------------------------------
   		! get the y-averaged density very near x=Lx/2
   		!----------------------------------------------
   		if( proc_row(YBLOCK,myid) == p1/2 ) then    ! all other procs just keep their zeros
    		i=1
    		do k=1,locnz
     			tmp(k) = SUM(s1(:,i,k))/locny       ! y avg density in [kg/m3] near Lx/2
    		enddo
   		endif
   		!------------------------------------------------------
   		! share values with all other procs with same proc_col
   		!------------------------------------------------------
   		do j=0,p2-1   
    		if( proc_col(YBLOCK,myid) == j ) then
     			call MPI_ALLREDUCE(tmp,rho_bar,locnz,MPI_DOUBLE_PRECISION,MPI_SUM,rows_comm(j),ierr)
    		endif    
   		enddo
   
   		deallocate(tmp)
		first_entry = .FALSE.
	endif
  
  
  
	!----------------------------------------
	!   average u in x within myid's YBLOCK
	!----------------------------------------
	id = 6  ! scratch
	X_means(:,:,id) = 0.d0
	do k=1,locnz
		do j=1,locny   
			do i=1,locnx
    			X_means(j,k,id) = X_means(j,k,id) + u(j,i,k)
    		enddo
   		enddo
  	enddo
  	X_means(:,:,id) = X_means(:,:,id) / dfloat(locnx)
  	!----------------------------------------
  	!   do the x avgs across processors
  	!----------------------------------------
  	id = 1   !  u
  	do j = 0,p2-1
  		if( proc_col(YBLOCK,myid) == j ) then
   			call MPI_ALLREDUCE(X_means(1,1,6),X_means(1,1,id),count,MPI_DOUBLE_PRECISION,MPI_SUM,rows_comm(j),ierr)
    		X_means(:,:,id) = X_means(:,:,id) / dfloat(p1)
   		endif
  	enddo

  
  
  	!----------------------------------------
  	!   average v in x within myid's YBLOCK
  	!----------------------------------------
  	id = 6  ! scratch
  	X_means(:,:,id) = 0.d0
  	do k=1,locnz
   		do j=1,locny   
    		do i=1,locnx
     			X_means(j,k,id) = X_means(j,k,id) + v(j,i,k)
    		enddo
   		enddo
  	enddo
  	X_means(:,:,id) = X_means(:,:,id) / dfloat(locnx)
  	!----------------------------------------
  	!   do the x avgs across processors
  	!----------------------------------------
  	id = 2   !  v
  	do j = 0,p2-1
  		if( proc_col(YBLOCK,myid) == j ) then
    		call MPI_ALLREDUCE(X_means(1,1,6),X_means(1,1,id),count,MPI_DOUBLE_PRECISION,MPI_SUM,rows_comm(j),ierr)
    		X_means(:,:,id) = X_means(:,:,id) / dfloat(p1)
   		endif
  	enddo
  
  	!----------------------------------------
  	!   average w in x within myid's YBLOCK
  	!----------------------------------------
  	id = 6  ! scratch
  	X_means(:,:,id) = 0.d0
  	do k=1,locnz
  		do j=1,locny   
    		do i=1,locnx
     			X_means(j,k,id) = X_means(j,k,id) + w(j,i,k)
    		enddo
   		enddo
  	enddo
  	X_means(:,:,id) = X_means(:,:,id) / dfloat(locnx)
  	!----------------------------------------
  	!   do the x avgs across processors
  	!----------------------------------------
  	id = 3   !  w
  	do j = 0,p2-1
   		if( proc_col(YBLOCK,myid) == j ) then
    		call MPI_ALLREDUCE(X_means(1,1,6),X_means(1,1,id),count,MPI_DOUBLE_PRECISION,MPI_SUM,rows_comm(j),ierr)
    		X_means(:,:,id) = X_means(:,:,id) / dfloat(p1)
   		endif
  	enddo
  
  	!----------------------------------------
  	!   average s1 in x within myid's YBLOCK
  	!----------------------------------------
  	id = 6  ! scratch
  	X_means(:,:,id) = 0.d0
  	do k=1,locnz
  		do j=1,locny   
    		do i=1,locnx
     			X_means(j,k,id) = X_means(j,k,id) + s1(j,i,k)
    		enddo
   		enddo
  	enddo
  	X_means(:,:,id) = X_means(:,:,id) / dfloat(locnx)
  	!----------------------------------------
  	!   do the x avgs across processors
  	!----------------------------------------
  	id = 4   !  s1
  	do j = 0,p2-1
		if( proc_col(YBLOCK,myid) == j ) then
    		call MPI_ALLREDUCE(X_means(1,1,6),X_means(1,1,id),count,MPI_DOUBLE_PRECISION,MPI_SUM,rows_comm(j),ierr)
    		X_means(:,:,id) = X_means(:,:,id) / dfloat(p1)
   		endif
  	enddo
  
  	!----------------------------------------
  	!   average s2 in x within myid's YBLOCK
  	!----------------------------------------
  	if( do_second_scalar ) then
   		id = 6  ! scratch
   		X_means(:,:,id) = 0.d0
   		do k=1,locnz
    		do j=1,locny   
     			do i=1,locnx
      				X_means(j,k,id) = X_means(j,k,id) + s2(j,i,k)
     			enddo
    		enddo
   		enddo
   		X_means(:,:,id) = X_means(:,:,id) / dfloat(locnx)
   		!----------------------------------------
   		!   do the x avgs across processors
   		!----------------------------------------
   		id = 5   !  s2
   		do j = 0,p2-1
    		if( proc_col(YBLOCK,myid) == j ) then
     			call MPI_ALLREDUCE(X_means(1,1,6),X_means(1,1,id),count,MPI_DOUBLE_PRECISION,MPI_SUM,rows_comm(j),ierr)
     			X_means(:,:,id) = X_means(:,:,id) / dfloat(p1)
    		endif
   		enddo
 	 endif
  
end subroutine user_analysis
