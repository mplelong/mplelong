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
	integer           :: id,endpt,i,j,k,idim
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
  	call mpi_bcast(y_periodic,1,MPI_LOGICAL,0,comm,ierr)
  	
  	if(myid==0) read(1,*) z_periodic
  	call mpi_bcast(z_periodic,1,MPI_LOGICAL,0,comm,ierr)
  
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
  		if(myid==0) read(1,*) p(3)
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
 	write_normal_derivs(:) = .FALSE.    ! only switch on as appropriate
 
 
 	do i=1,num_file_sets
  		if(myid==0) read(1,*) filename_root(i)
  		call mpi_bcast(filename_root(i),80,MPI_CHARACTER,0,comm,ierr)
  
  		if(myid==0) read(1,*) mode(i)
  		call mpi_bcast(mode(i),80,MPI_CHARACTER,0,comm,ierr)
  
  		if(myid==0) read(1,*) nsteps(i)
  		call mpi_bcast(nsteps(i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) ilocs(1,i)
  		if( ilocs(1,i) < 0 ) ilocs(1,i) = 1    ! negative values trigger limiting value 
  		call mpi_bcast(ilocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) ilocs(2,i)
  		if( ilocs(2,i) < 0 ) ilocs(2,i) = nx   ! negative values trigger limiting value 
  		call mpi_bcast(ilocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) ilocs(3,i)
  		call mpi_bcast(ilocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) jlocs(1,i)
  		if( jlocs(1,i) < 0 ) jlocs(1,i) = 1    ! negative values trigger limiting value
  		call mpi_bcast(jlocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) jlocs(2,i)
  		if( jlocs(2,i) < 0 ) jlocs(2,i) = ny   ! negative values trigger limiting value
  		call mpi_bcast(jlocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) jlocs(3,i)
  		call mpi_bcast(jlocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) klocs(1,i)
  		if( klocs(1,i) < 0 ) klocs(1,i) = 1    ! negative values trigger limiting value
  		call mpi_bcast(klocs(1,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) klocs(2,i)
  		if( klocs(2,i) < 0 ) klocs(2,i) = nz   ! negative values trigger limiting value
  		call mpi_bcast(klocs(2,i),1,MPI_INTEGER,0,comm,ierr)
  
  		if(myid==0) read(1,*) klocs(3,i)
  		call mpi_bcast(klocs(3,i),1,MPI_INTEGER,0,comm,ierr)
  
  		do j=1,5 !! u,v,w,s1,s2
   			if(myid==0) read(1,*) variable_key(j,i)
   			call mpi_bcast(variable_key(j,i),1,MPI_INTEGER,0,comm,ierr)   ! was LOGIVCAL  why did this even work???
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
 	
 	!   read logical variable to see if child grid output needs to be configured
 	if(myid==0) read(1,*) do_child_grid
   	call mpi_bcast(do_child_grid,1,MPI_LOGICAL,0,comm,ierr) 
   	
   	if( do_child_grid ) then
   	
   		!--------------------------------------------------------
   		!  read the extra child grid parameters
   		!--------------------------------------------------------
   		if(myid==0) read(1,*) i0_child
  		call mpi_bcast(i0_child,1,MPI_INTEGER,0,comm,ierr)
  		if(myid==0) read(1,*) i1_child
  		call mpi_bcast(i1_child,1,MPI_INTEGER,0,comm,ierr)
  		
  		if(myid==0) read(1,*) j0_child
  		call mpi_bcast(j0_child,1,MPI_INTEGER,0,comm,ierr)
  		if(myid==0) read(1,*) j1_child
  		call mpi_bcast(j1_child,1,MPI_INTEGER,0,comm,ierr)
  		
  		if(myid==0) read(1,*) k0_child
  		call mpi_bcast(k0_child,1,MPI_INTEGER,0,comm,ierr)
  		if(myid==0) read(1,*) k1_child
  		call mpi_bcast(k1_child,1,MPI_INTEGER,0,comm,ierr)
  		
  		if(myid==0) read(1,*) inc_child_2d
  		call mpi_bcast(inc_child_2d,1,MPI_INTEGER,0,comm,ierr)
  		
  		if(myid==0) read(1,*) inc_child_3d
  		call mpi_bcast(inc_child_3d,1,MPI_INTEGER,0,comm,ierr)
  		
  		if(myid==0) read(1,*) t_start_child
		call mpi_bcast(t_start_child,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
		
		if(myid==0) read(1,*) write_s2_child_grid
		call mpi_bcast(write_s2_child_grid,1,MPI_LOGICAL,0,comm,ierr)
		
		if(myid==0) read(1,*) write_s1_bar_child_grid
		call mpi_bcast(write_s1_bar_child_grid,1,MPI_LOGICAL,0,comm,ierr)
		
		if(myid==0) read(1,*) write_s2_bar_child_grid
		call mpi_bcast(write_s2_bar_child_grid,1,MPI_LOGICAL,0,comm,ierr)
		!--------------------------------------------------------
   		!  done reading the extra child grid parameters
   		!--------------------------------------------------------
		
		
		!----------------------------------------------------------
		! Add 7 to num_file_sets for E/W, S/N, B/T and XYZ_child
		!----------------------------------------------------------
   		j = num_file_sets    
   		num_file_sets = num_file_sets + 7
   		
   		!-------------------------------------------------------------
   		! indicate that these files need normal derivs to be written
   		!-------------------------------------------------------------
   		write_normal_derivs(i+1:i+7) = .TRUE. 
   		
   		!--------------------
   		! set the file names
   		!--------------------
   		filename_root(j+1) = 'east'
   		filename_root(j+2) = 'west'
   		filename_root(j+3) = 'south'
   		filename_root(j+4) = 'north'
   		filename_root(j+5) = 'bottom'
   		filename_root(j+6) = 'top'
   		filename_root(j+7) = 'XYZ_child'
   		
   		!-----------------------------------------------
   		! all processors set the write parameters
   		! afterwards, these should be look like all
   		! the other filesets except for the additional
   		! output of the normal derivatives
   		!-----------------------------------------------
   		if(myid==0) open(11,file='output/debug_data/child_grid_data')
   		do j=num_file_sets-7+1,num_file_sets   ! j is fileset id
   			
   			variable_key(6,j) = 0                ! suppress output of div_u*
			variable_key(9:11,j) = 0             ! suppress  output of ustar,vstar and wstar
		 			
  			if(j<num_file_sets) then
  				nsteps(j) = inc_child_2d
  				mode(j) = 'append'
  			else
  				nsteps(j) = inc_child_3d
  				mode(j) = 'new'
  			endif
  			
  			ilocs(1,j) = i0_child
  			ilocs(2,j) = i1_child
  			ilocs(3,j) = 1
  			jlocs(1,j) = j0_child
  			jlocs(2,j) = j1_child
  			jlocs(3,j) = 1
  			klocs(1,j) = k0_child
  			klocs(2,j) = k1_child
  			klocs(3,j) = 1
  			
  			k = num_file_sets - 7
  			if( j==k+1 )  ilocs(2,j) = i0_child   !  E
  			if( j==k+2 )  ilocs(1,j) = i1_child   !  W
  			if( j==k+3 )  jlocs(2,j) = j0_child   !  S
  			if( j==k+4 )  jlocs(1,j) = j1_child   !  N
  			if( j==k+5 )  klocs(2,j) = k0_child   !  B
  			if( j==k+6 )  klocs(1,j) = k1_child   !  T
  			
  			
  			variable_key(1:4,j) = 1            	! always write u,v,w,s1  			  			
  			if( write_s2_child_grid ) then
  				variable_key(5,j) = 1     		! write s2
  			else
  				variable_key(5,j) = 0
  			endif
  			
  			if(myid==0) then
  				write(11,*) ' '
  				write(11,*) '         fileset , root: ',j,filename_root(j)
  				write(11,*) '                  mode : ', mode(j)
  				write(11,*) '                    inc: ', nsteps(j)
  				write(11,*) '          t_start_child: ', t_start_child
  				write(11,*) '                  ilocs: ', ilocs(:,j)
  				write(11,*) '                  jlocs: ', jlocs(:,j)
  				write(11,*) '                  klocs: ', klocs(:,j) 				
  				write(11,*) '     variable_key(:,j) : ', variable_key(1:5,j)  				
  				write(11,*) '    write_normal_derivs: ', write_normal_derivs(j)
  				write(11,*) '    write_s2_child_grid: ', write_s2_child_grid
  				write(11,*) 'write_s1_bar_child_grid: ', write_s1_bar_child_grid
  				write(11,*) 'write_s2_bar_child_grid: ', write_s2_bar_child_grid
  				write(11,*) ' '
  				write(11,*) ' '
  			endif
  			  			
  		enddo
  		if(myid==0) close(11)
  		 		
   	endif ! end do_child_grid block
   	
 		
 	if(myid==0) write(0,*) ' ................      ',trim(file2),' read'
 	if(myid==0) close(1)


	call mpi_barrier(comm,ierr)
 return
end subroutine ReadUserData
