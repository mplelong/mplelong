!-------------------------------------------------------------!
!  flow_solve                                   4/2021        !
!                                                             !
!    MPI parallel version, f90                                !
!    Two-dimensional data and processor decomposition         !
!                                                             !
!    Spatial differentiation via Bernoulli-cos method with    !
!    no particular symmetries assumed or imposed at bdries    !
!    Boundary condition data supplied by user                 ! 
!       (i.e. from a parent run)                              !
!                                                             !
!    1 or 2 scalars, configured as combinations of density,   !
!    temperature, salinity or passive tracers.                !
!                                                             !
!    2nd, 3rd or 4th order Adams Bashforth time stepping      !
!    for advection, buoyancy, rotation, and forcing.          !
!                                                             !
!    Exact integration of arbitrary-even-order diffusion      !    
!    operators in wavenumber space                            !
!                                                             !
!    Hooks for user-defined routines to implement             !
!    xtra terms on rhs's, i.e. forcing, sponge layers etc     !
!                                                             !
!    External software for flow_solve                         !
!     MPI               (message passing)                     !
!     BLAS              (basic linear algebra, lapack)        !
!     FFTW3             (sequential spectral transforms)      !
!     NETCDF            (output data format)                  !
!                                                             !
!    Useful software for data processing etc                  !
!     NCO               (netcdf operators)                    !
!     python            (mpi4py,netCDF4,scipy,matplotlib)     !
!     perl              (w/ nco for manipulating .nc files)   !
!     ncview            (quick viewing of .nc files)          !
!                                                             !
!    Kraig Winters                                            !
!    Scripps Institution of Oceanography                      !
!    University of California, San Diego                      !
!-------------------------------------------------------------!  
program flow_solve
use etc
call preliminary_tasks
 do while (integrate);
!  call user_analysis;
   call explicit_rhs;STOP "B"
    call apply_forcing;STOP "C"
     call explicit_step;STOP "D"
      call fill_boundary_arrays;STOP "E"
       call pressure_projection;STOP "F"
        call apply_bcs;STOP "G"
         call diffuse        
          call write_results
           call toggle
            enddo
end program flow_solve


