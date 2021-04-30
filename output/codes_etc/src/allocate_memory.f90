subroutine AllocateIndependentVariables
	use etc,                   only: logfile
	use mpi_params,            only: myid,comm,ierr 
	use independent_variables, only: x,y,z,nx,ny,nz   ! global coord arrays
	implicit none
	allocate( x(nx), y(ny), z(nz) )
end subroutine AllocateIndependentVariables


subroutine AllocateDependentVariables
	use etc,                   only: logfile
	use mpi_params,            only: myid,comm,ierr 
	use decomposition_params 
	use dependent_variables
	use independent_variables, only: nz
	use methods_params
	implicit none 
 
	integer                       :: npvs,i,n
	real(kind=8)                  :: words_in_yblock
	real(kind=8),parameter        :: dzero=0.d0

 
	if(myid==0) then
		write(0,*) ' ................'
		write(0,*) ' ................     hello world from AllocateDependentVariables'
		open(1,file=logfile,position='append') 
		write(1,*) '  '
		write(1,*) '  '
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '                  AllocateDependentVariables Report:'
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '  '
	endif
 
 	!! Allocate space for dependent variables in YBLOCK decomposition
 
 	words_in_yblock = array_size(IDIM,YBLOCK,myid)           &
                     *array_size(JDIM,YBLOCK,myid)           &
                     *array_size(KDIM,YBLOCK,myid)
  
                    
 	allocate(  u( array_size(IDIM,YBLOCK,myid),              &
                  array_size(JDIM,YBLOCK,myid),              &
                  array_size(KDIM,YBLOCK,myid)  ) )
 
 	allocate(  v( array_size(IDIM,YBLOCK,myid),              &
                  array_size(JDIM,YBLOCK,myid),              &
                  array_size(KDIM,YBLOCK,myid)  ) )
 
 	allocate(  w( array_size(IDIM,YBLOCK,myid),              &
                  array_size(JDIM,YBLOCK,myid),              &
                  array_size(KDIM,YBLOCK,myid)  ) )
 
 	allocate( s1( array_size(IDIM,YBLOCK,myid),              &
                  array_size(JDIM,YBLOCK,myid),              &
                  array_size(KDIM,YBLOCK,myid)  ) )
 	allocate ( s1_bar(nz,3) )
 	n = words_in_yblock
 	call dinit(n,dzero,s1)
 	n = 3*nz
 	call dinit(n,dzero,s1_bar)
 
 	if( do_second_scalar ) then
  		allocate( s2( array_size(IDIM,YBLOCK,myid),              &
                      array_size(JDIM,YBLOCK,myid),              &
                     array_size(KDIM,YBLOCK,myid)  ) )
  		allocate ( s2_bar(nz,3) )
  		n = words_in_yblock
  		call dinit(n,dzero,s2)
  		n = 3*nz
  		call dinit(n,dzero,s1_bar)
  		npvs=5
 	else
  		npvs=4
 	endif
 
 	n = words_in_yblock 
 	call dinit(n,dzero,u)
 	call dinit(n,dzero,v)
 	call dinit(n,dzero,w)
 
 
 	if(myid==0) then
  		do i=0,1
			write(i,*) ' ................      local size, u array  in Mb                   ',words_in_yblock*8./(1024.d0*1024.)
			write(i,*) ' ................      local size, all primitive variables  in Mb   ',npvs*words_in_yblock*8./(1024.*1024.)
			write(i,*) ' ................'
		enddo
		write(1,*) ' -----> AllocateDependentVariables routine exiting normally  <---------- '
  		close(1)
	endif

 	call mpi_barrier(comm,ierr)
 	return
end subroutine AllocateDependentVariables



subroutine AllocateIntermediateVariables
	use etc,                   only: logfile
	use mpi_params,            only: myid,comm,ierr,numprocs
	use independent_variables, only: nx,ny
	use decomposition_params 
	use intermediate_variables
	use methods_params
	implicit none 
	integer                       :: narrays(3),npvs,i,n,ksize
	real(kind=8)                  :: words_in_xblock
	real(kind=8)                  :: words_in_yblock
	real(kind=8)                  :: words_in_zblock
	real(kind=8),parameter        :: dzero=0.d0
 
 
	if(myid==0) then
		write(0,*) ' ................'
		write(0,*) ' ................     hello world from AllocateIntermediateVariables'
		open(1,file=logfile,position='append') 
		write(1,*) '  '
		write(1,*) '  '
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '                  AllocateIntermediateVariables Report:'
		write(1,*) ' =========================================================== '
		write(1,*) ' =========================================================== '
		write(1,*) '  '
	endif
 
 	!-----------------------------------------------------------------
 	!! Allocate space for work arrays in each BLOCK decomposition
 	!-----------------------------------------------------------------
 	narrays(:)=(/3,6,3/)   ! needed extra for div grad approach to Laplacian diffusion
 
 	if( do_second_scalar ) then
		npvs=5
	else
		npvs=4
	endif

 
 	words_in_xblock = array_size(IDIM,XBLOCK,myid)              &
                     *array_size(JDIM,XBLOCK,myid)              &
                     *array_size(KDIM,XBLOCK,myid)
 
 	words_in_yblock = array_size(IDIM,YBLOCK,myid)              &
                     *array_size(JDIM,YBLOCK,myid)              &
                     *array_size(KDIM,YBLOCK,myid)
                  
 	words_in_zblock = array_size(IDIM,ZBLOCK,myid)              &
                     *array_size(JDIM,ZBLOCK,myid)              &
                     *array_size(KDIM,ZBLOCK,myid)
                  
 	allocate(  tmpX( array_size(IDIM,XBLOCK,myid),              &
                     array_size(JDIM,XBLOCK,myid),              &
                     array_size(KDIM,XBLOCK,myid),narrays(1)  ))
  
 	allocate(  tmpY( array_size(IDIM,YBLOCK,myid),              &
                     array_size(JDIM,YBLOCK,myid),              &
                     array_size(KDIM,YBLOCK,myid),narrays(2) ))
                  
 	allocate(  tmpZ( array_size(IDIM,ZBLOCK,myid),              &
                     array_size(JDIM,ZBLOCK,myid),              &
                     array_size(KDIM,ZBLOCK,myid),narrays(3) ))
 
  
  	n = words_in_xblock*narrays(1)
  	call dinit(n,dzero,tmpX)
  	n = words_in_yblock*narrays(2)
  	call dinit(n,dzero,tmpY)
  	n = words_in_zblock*narrays(3)
  	call dinit(n,dzero,tmpZ)

  
	!-----------------------------------------------------------------------
	!  Allocate space for rhs arrays in YBLOCK format
	!-----------------------------------------------------------------------
  
	allocate(  explicit_rhs( array_size(IDIM,YBLOCK,myid),              &
                             array_size(JDIM,YBLOCK,myid),              &
                             array_size(KDIM,YBLOCK,myid),              &
                             npvs,AB_ORDER))  

  	n = words_in_yblock*npvs*(AB_ORDER)
  	call dinit(n,dzero,explicit_rhs)

                   
	!-----------------------------------------------------------------------
	!  Allocate space for phi and the rhs for elliptic eqn for phi
	!-----------------------------------------------------------------------
  
	allocate( phi( array_size(IDIM,YBLOCK,myid),              &
                   array_size(JDIM,YBLOCK,myid),              &
                   array_size(KDIM,YBLOCK,myid) ) )
  
	allocate(div_u( array_size(IDIM,YBLOCK,myid),              &
                    array_size(JDIM,YBLOCK,myid),              &
                    array_size(KDIM,YBLOCK,myid) ) )
                    
    allocate(ustar( array_size(IDIM,YBLOCK,myid),              &
                    array_size(JDIM,YBLOCK,myid),              &
                    array_size(KDIM,YBLOCK,myid) ) )
    
    allocate(vstar( array_size(IDIM,YBLOCK,myid),              &
                    array_size(JDIM,YBLOCK,myid),              &
                    array_size(KDIM,YBLOCK,myid) ) )
	
	allocate(wstar( array_size(IDIM,YBLOCK,myid),              &
                    array_size(JDIM,YBLOCK,myid),              &
                    array_size(KDIM,YBLOCK,myid) ) )
                    
  	n = words_in_yblock
  	call dinit(n,dzero,phi)
  	call dinit(n,dzero,div_u)
  	call dinit(n,dzero,ustar)
  	call dinit(n,dzero,vstar)
  	call dinit(n,dzero,wstar)
  
  	!-----------------------------------------------------------------------
  	!  Allocate space for reading in YZ data for initializing 3d runs
  	!-----------------------------------------------------------------------
 	allocate(  tmpYZ( array_size(IDIM,YBLOCK,myid),              &
			                                     1,              &
					  array_size(KDIM,YBLOCK,myid)  ) )
  	n = array_size(IDIM,YBLOCK,myid)*array_size(KDIM,YBLOCK,myid)
  	call dinit(n,dzero,tmpYZ)
  

 	if(myid==0) then
  		do i=0,1
			write(i,*) ' ................      local size, tmpX array(s) in Mb             ',words_in_xblock*8./(1024.*1024.)*narrays(1)
   			write(i,*) ' ................      local size, tmpY array(s) in Mb             ',words_in_yblock*8./(1024.*1024.)*narrays(2)
   			write(i,*) ' ................      local size, tmpZ array(s) in Mb             ',words_in_zblock*8./(1024.*1024.)*narrays(3)
   			write(i,*) ' ................      local size, explicit_rhs array(s) Mb        ',words_in_yblock*8./(1024.*1024.)*npvs*(AB_ORDER)
   			write(i,*) ' ................      local size, phi array in Mb                 ',words_in_yblock*8./(1024.*1024.)
   			write(i,*) ' ................      local size, div_u array in Mb               ',words_in_yblock*8./(1024.*1024.)
   			write(i,*) ' ................      local size, ustar array in Mb               ',words_in_yblock*8./(1024.*1024.)
   			write(i,*) ' ................      local size, vstar array in Mb               ',words_in_yblock*8./(1024.*1024.)
   			write(i,*) ' ................      local size, wstar array in Mb               ',words_in_yblock*8./(1024.*1024.)
   			write(i,*) ' ................'
  		enddo
  		write(1,*) ' -----> AllocateIntermediateVariables routine exiting normally  <---------- '
  		close(1)
 	endif

 	call mpi_barrier(comm,ierr)
 	return
end subroutine AllocateIntermediateVariables


