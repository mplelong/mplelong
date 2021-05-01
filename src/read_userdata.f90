subroutine ReadUserData
	use etc
	use dimensional_scales
	use dimensionless_params
	use mpi_params
	use methods_params
	use decomposition_params,   only: np,p1,p2 
	use io_params  
	use dependent_variables
	use independent_variables
	use intermediate_variables
	implicit none  

	include 'mpif.h' 
	integer           :: id,endpt,i,j,idim
	character(len=80) :: file1,file2,file3,file4,file5

	if(myid==0) write(0,*) ' ................'
	if(myid==0) write(0,*) ' ................      hello world from read_user_data'
 

	file1='input/problem_params'
	file2='input/io_params'

  
	if(myid==0) open(1,file=file1,position='rewind') 
  
	if(myid==0) read(1,*) runlabel
	call mpi_bcast(runlabel,80,MPI_CHARACTER,0,comm,ierr)
	
	if(myid==0) read(1,*) restart
	call mpi_bcast(restart,1,MPI_LOGICAL,0,comm,ierr)
          
	if(myid==0) read(1,*) do_second_scalar
	call mpi_bcast(do_second_scalar,1,MPI_LOGICAL,0,comm,ierr)
    
    if(myid==0) read(1,*) AB_order
	call mpi_bcast(AB_order,1,MPI_INTEGER,0,comm,ierr)
	      
	if(myid==0) read(1,*) p1
	call mpi_bcast(p1,1,MPI_INTEGER,0,comm,ierr)
  
	if(myid==0) read(1,*) p2
	call mpi_bcast(p2,1,MPI_INTEGER,0,comm,ierr)
  
	if(myid==0) read(1,*) nx
	call mpi_bcast(nx,1,MPI_INTEGER,0,comm,ierr)
  
	if(myid==0) read(1,*) ny
	call mpi_bcast(ny,1,MPI_INTEGER,0,comm,ierr)
  
	if(myid==0) read(1,*) nz
	call mpi_bcast(nz,1,MPI_INTEGER,0,comm,ierr)
  
	if(myid==0) read(1,*) dt
	call mpi_bcast(dt,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
	if(myid==0) read(1,*) t0
	call mpi_bcast(t0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
	if(myid==0) read(1,*) tf
	call mpi_bcast(tf,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
	if(myid==0) read(1,*) Lx
	call mpi_bcast(Lx,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
	if(myid==0) read(1,*) Ly
	call mpi_bcast(Ly,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
	if(myid==0) read(1,*) Lz
	call mpi_bcast(Lz,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
	
	if(myid==0) read(1,*) x_periodic
  	call mpi_bcast(x_periodic,1,MPI_LOGICAL,0,comm,ierr)
  	
  	if(myid==0) read(1,*) y_periodic
  	call mpi_bcast(x_periodic,1,MPI_LOGICAL,0,comm,ierr)
  	
  	if(myid==0) read(1,*) z_periodic
  	call mpi_bcast(x_periodic,1,MPI_LOGICAL,0,comm,ierr)
  
  	if(myid==0) read(1,*) scalar_kind(1)
  	call mpi_bcast(scalar_kind(1),1,MPI_CHARACTER,0,comm,ierr)
  
  	if(myid==0) read(1,*) scalar_kind(2)
  	call mpi_bcast(scalar_kind(2),1,MPI_CHARACTER,0,comm,ierr)
  
  	if(myid==0) read(1,*) do_forcing
  	call mpi_bcast(do_forcing,1,MPI_LOGICAL,0,comm,ierr)
        
	if(myid==0) read(1,*) rho_0
	call mpi_bcast(rho_0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
	rho0 = rho_0
	
	if(myid==0) read(1,*) g
	call mpi_bcast(g,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
	
	if(myid==0) read(1,*) f0
	call mpi_bcast(f0,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
	coriolis = f0
  
	if(myid==0) read(1,*) nu
	call mpi_bcast(nu,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  	if(myid==0) read(1,*) kappa(1)
  	call mpi_bcast(kappa(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  	if(myid==0) read(1,*) kappa(2)
  	call mpi_bcast(kappa(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  
  	if(myid==0) read(1,*) high_order_operators
  	call mpi_bcast(high_order_operators,1,MPI_LOGICAL,0,comm,ierr)
  	
  	!  if necessary, keep reading to get the params for high-order diffusion
  	if( high_order_operators ) then
  		if(myid==0) read(1,*) p(1)
  		call mpi_bcast(p(1),1,MPI_INTEGER,0,comm,ierr)
  		if(myid==0) read(1,*) p(2)
  		call mpi_bcast(p(2),1,MPI_INTEGER,0,comm,ierr)
  		if(myid==0) read(1,*) p(2)
  		call mpi_bcast(p(3),1,MPI_INTEGER,0,comm,ierr)
  		
  		if(myid==0) read(1,*) T_diff(1)
  		call mpi_bcast(T_diff(1),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  		if(myid==0) read(1,*) T_diff(2)
  		call mpi_bcast(T_diff(2),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  		if(myid==0) read(1,*) T_diff(3)
  		call mpi_bcast(T_diff(3),1,MPI_DOUBLE_PRECISION,0,comm,ierr)
  	endif
  
  
	if(myid==0) close(1)
 	if(np .ne. p1*p2)       stop ' processor grid error, p1*p2 NE np '     
 	if(numprocs .ne. np)    stop ' processor grid error, numprocs NE np '    
 	if(myid==0) write(0,*) ' ................      ',trim(file1),' read'

 
	if(myid==0) open(1,file=file2,position='rewind')
	if(myid==0) read(1,*) num_file_sets
 	call mpi_bcast(num_file_sets,1,MPI_INTEGER,0,comm,ierr)
 
 
 	do i=1,num_file_sets
  		if(myid==0) read(1,*) filename_root(i)
  		call mpi_bcast(filename_root(i),80,MPI_CHARACTER,0,comm,ierr)
  
  		if(myid==0) read(1,*) mode(i)
  		call mpi_bcast(mode(i),80,MPI_CHARACTER,0,comm,ierr)
  
  		if(myid==0) read(1,*) nsteps(i)
  		call mpi_bcast(nsteps(i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) ilocs(1,i)
  		call mpi_bcast(ilocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) ilocs(2,i)
  		call mpi_bcast(ilocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) ilocs(3,i)
  		call mpi_bcast(ilocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) jlocs(1,i)
  		call mpi_bcast(jlocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) jlocs(2,i)
  		call mpi_bcast(jlocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) jlocs(3,i)
  		call mpi_bcast(jlocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) klocs(1,i)
  		call mpi_bcast(klocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) klocs(2,i)
  		call mpi_bcast(klocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) klocs(3,i)
  		call mpi_bcast(klocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  		do j=1,5 !! u,v,w,s1,s2
   			if(myid==0) read(1,*) variable_key(j,i)
   			call mpi_bcast(variable_key(j,i),1,MPI_LOGICAL,0,comm,ierr)
  		enddo
  
  		if( variable_key(4,i) == 1 ) write_s1_bar(i)=.FALSE.
  		if( variable_key(4,i) == 2 ) write_s1_bar(i)=.TRUE.
  
  		if( variable_key(5,i) == 1 ) write_s2_bar(i)=.FALSE.
  		if( variable_key(5,i) == 2 ) write_s2_bar(i)=.TRUE.
  
  		if(myid==0) write(0,*) ' ................      output fileset number ',i,' read'
  
  		!! set these to 1 during testing; generally 0 to save disk space
  		variable_key(6,i)=0   ! divustar
  
  		!! set permanently to 1
  		variable_key(7,i)=1   ! phi
  
  		!! not pd though
  		variable_key(8,i)=0   ! pd
  
  		!! quick sanity check ....
  		if(ilocs(2,i) > nx ) stop 'io_params data > nx'
  		if(jlocs(2,i) > ny ) stop 'io_params data > ny'
  		if(klocs(2,i) > nz ) stop 'io_params data > nz'
 	enddo 
 	if(myid==0) write(0,*) ' ................      ',trim(file2),' read'
 	if(myid==0) close(1)


	call mpi_barrier(comm,ierr)
	return
end subroutine ReadUserData
