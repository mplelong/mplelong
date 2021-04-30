SHELL=/bin/sh
VENDOR=intel

#--------------------------------------------------------------------
#    Root directories for external software... generally
#    below these directories will be include, lib, bin
#--------------------------------------------------------------------
MPI_ROOT= 
FFTW3_ROOT=

#--------------------------------------------------------------------
#    Run and output base directory locations
#--------------------------------------------------------------------
BASE=${PWD}

#--------------------------------------------------------------------
#    Compiler wrappers for MPI
#    intel mpi enabled by default
#--------------------------------------------------------------------
F77 = mpif77 
F90 = mpif90  
FLINKER = ${F90}
 
#--------------------------------------------------------------------
#    Vendor specific compilation parameters
#--------------------------------------------------------------------

ifeq ($(VENDOR), portland)
    NETCDF_ROOT=/opt/netcdf/4.3.2/pgi/mvapich2_ib
	LDFLAGS = 
	FOPTFLAGS   = -Mpreprocess -Mextend -g -fast -Mipa=fast,inline #-Mbounds -Mchkptr    #-Ktrap=fp
	F77OPTFLAGS = -Mpreprocess -Mextend -g -fast -Mipa=fast,inline #-Mbounds             #-Ktrap=fp
	BLAS_LIB	= -L/home/diag/opt/atlas/pgi/lib -lblas
endif
ifeq ($(VENDOR), intel)
    NETCDF_ROOT=/opt/netcdf/4.3.2/intel/mvapich2_ib
	LDFLAGS = 
	FOPTFLAGS   = -xHost -132 -O3 -static -ip -no-prec-div -static-intel #-check bounds -traceback
	F77OPTFLAGS = -xHost -132 -O3 -static -ip -no-prec-div -static-intel #-check bounds -traceback
	BLAS_LIB	= -L/opt/atlas/intel/lib -lf77blas -latlas
endif
# -132 -O1 -ip -opt-prefetch -assume buffered_io #-check bounds -traceback 

XTRA_LIBS =			#  in case you need something special on your machine


#--------------------------------------------------------------------
#  include and link flags...
#--------------------------------------------------------------------  
NETCDF_INC=-I${NETCDF_ROOT}/include
NETCDF_LIB=-L${NETCDF_ROOT}/lib -lnetcdf -lnetcdff

FFTW_INC=-I${FFTW3_ROOT}/include
FFTW_LIB=-L${FFTW3_ROOT}/lib -lfftw3

ALL_EXTERNAL_LIBS =	$(NETCDF_LIB) $(FFTW_LIB) $(BLAS_LIB) $(XTRA_LIBS) 		
ALL_INCLUDES =  $(NETCDF_INC) $(FFTW_INC)	

 
showconfig:
	@echo
	-@echo  "--------------------------------------------------------------------------------------" 
	-@echo  "host:                 " `hostname -s`
	-@echo  "OS:                   " `uname -sr`
	-@echo  "user:                 " `whoami`
	-@echo  "compile date:         " `date -u +"%Y-%m-%dT%H:%M:%SZ"`
	-@echo  "SHELL:                " ${SHELL}
	-@echo  "VENDOR:               " ${VENDOR}
	-@echo  "LD_LIBRARY_PATH:      " ${LD_LIBRARY_PATH}
	-@echo  "DYLD_LIBRARY_PATH:    " ${DYLD_LIBRARY_PATH}
	-@echo  "BASE:                 " ${BASE}
	-@echo  "OUTPUT_ROOT:          " ${BASE}/output
	-@echo  "F90:                  " ${F90}
	-@echo  "F77:                  " ${F77}
	-@echo  "FOPTFLAGS:            " ${FOPTFLAGS}
	-@echo  "F77OPTFLAGS:          " ${F77OPTFLAGS}
	-@echo  "FLINKER:              " ${FLINKER}
	-@echo  "NETCDF_INC:           " ${NETCDF_INC}
	-@echo  "NETCDF_LIB:           " ${NETCDF_LIB}
	-@echo  "FFTW_INC:             " ${FFTW_INC}
	-@echo  "FFTW_LIB:             " ${FFTW_LIB}
	-@echo  "BLAS_LIB:             " ${BLAS_LIB}
	-@echo  "ALL_EXTERNAL_LIBS:    " ${ALL_EXTERNAL_LIBS}
	-@echo  "ALL_INCLUDES:         " ${ALL_INCLUDES}
	-@echo  "XTRA_LIBS:            " ${XTRA_LIBS}
	-@echo  "--------------------------------------------------------------------------------------" 
