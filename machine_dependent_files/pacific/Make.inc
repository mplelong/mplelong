SHELL=/bin/sh
VENDOR=intel

#--------------------------------------------------------------------
#    Root directories for external software... 
#    (1) the env variables MPI_HOME and FFTW_ROOT are defined and
#        available via
#         module load intel
#         module load mpi/openmpi-3.1.0-intel
#         module load fftw
#    (2)  explicit netcdf path to local KW install, will pick up static libs
#         needed for compilation
#    (3)  at runtime (i.e. in the run.sh script)
#         I add the paths to the shared netcdf libraries needed by flow.x at runtime
#         export LD_LIBRARY_PATH= \
# $LD_LIBRARY_PATH:/home/kucla/KBW/kbw_installs/netcdf-c/lib:/home/kucla/KBW/kbw_installs/netcdf-fortran/lib
#         
#--------------------------------------------------------------------
MPI_ROOT= ${MPI_HOME}
FFTW3_ROOT=${FFTW_ROOT}
NETCDF_ROOT=${HOME}/KBW/kbw_installs/netcdf-fortran

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
ifeq ($(VENDOR), gcc)
	LDFLAGS = 
	FOPTFLAGS   = -g -O3 -static -cpp -ffree-form -ffree-line-length-none    #-fbounds-check -ffpe-trap=zero,overflow,underflow
	F77OPTFLAGS = -g -O3 -static -cpp -ffixed-form -ffixed-line-length-none  #-fbounds-check -ffpe-trap=zero,overflow,underflow
endif

ifeq ($(VENDOR), intel)
	LDFLAGS = 
	FOPTFLAGS   = -xHost -132 -O3 -ip -no-prec-div #-check bounds -traceback
	F77OPTFLAGS = -xHost -132 -O3 -ip -no-prec-div #-check bounds -traceback
endif

XTRA_LIBS = 	#  in case you need something special on your machine

#   some stray intel flags I've used in the past...
# -132 -O1 -ip -opt-prefetch -assume buffered_io  
# -static -static-intel




#--------------------------------------------------------------------
#  include and link flags...
# ${NETCDF_ROOT}/lib/libnetcdff.so.7 needed at runtime
# path added to LD_LIBRARY_PATH in ~/.bashrc and in run.sh
#--------------------------------------------------------------------  
NETCDF_INC= -I${NETCDF_ROOT}/include
NETCDF_LIB= -L${NETCDF_ROOT}/lib -lnetcdff 

FFTW_INC= -I${FFTW3_ROOT}/include
FFTW_LIB= -L${FFTW3_ROOT}/lib -lfftw3

BLAS_LIB= -L${BLAS_LIBRARY} -lf77blas -latlas

ALL_EXTERNAL_LIBS = $(FFTW_LIB) $(BLAS_LIB) $(NETCDF_LIB) $(XTRA_LIBS) 		
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
