SHELL=/bin/sh

#VENDOR=gcc_debug
VENDOR=gcc


#--------------------------------------------------------------------
#    Root directories for external software... generally
#    below these directories will be include, lib, bin
#
#        heron
#
#--------------------------------------------------------------------
MPI_ROOT = 
NETCDF_ROOT = /usr
FFTW3_ROOT = /usr

#--------------------------------------------------------------------
#    Run and output base directory locations
#--------------------------------------------------------------------
BASE=${PWD}

#--------------------------------------------------------------------
#    Compiler wrappers for MPI
#--------------------------------------------------------------------
F77 = mpif77
F90 = mpif90
FLINKER = ${F90}

#--------------------------------------------------------------------
#    Vendor specific compilation parameters
#    -fallow-argument-mismatch 
#--------------------------------------------------------------------

ifeq ($(VENDOR), gcc_debug)
	LDFLAGS =
	FOPTFLAGS   = -g -O3 -cpp -ffree-form -ffree-line-length-none -fallow-argument-mismatch -fbounds-check -ffpe-trap=zero,overflow,underflow -fbacktrace
	F77OPTFLAGS = -g -O3 -cpp -ffixed-form -ffixed-line-length-none -fallow-argument-mismatch -fbounds-check -ffpe-trap=zero,overflow,underflow -fbacktrace
	BLAS_LIB =  /usr/lib64/libblas.a
	LAPACK_LIB = /usr/lib64/liblapack.a
endif

ifeq ($(VENDOR), gcc)
	LDFLAGS =
	FOPTFLAGS   = -g -O3 -cpp -ffree-form -ffree-line-length-none -fallow-argument-mismatch 
	F77OPTFLAGS = -g -O3 -cpp -ffixed-form -ffixed-line-length-none -fallow-argument-mismatch 
	BLAS_LIB =  /usr/lib64/libblas.a
	LAPACK_LIB = /usr/lib64/liblapack.a
endif



XTRA_LIBS =						# in case you need something special on your machine
LD_LIBRARY_PATH =				# only search for explicitly indicated libraries
DYLD_LIBRARY_PATH =				# only search for explicitly indicated libraries


#--------------------------------------------------------------------
#  include and link flags...
#--------------------------------------------------------------------
NETCDF_INC = -I${NETCDF_ROOT}/include
NETCDF_LIB = -L${NETCDF_ROOT}/lib64 -lnetcdff -lnetcdf

FFTW_INC = -I${FFTW3_ROOT}/include
FFTW_LIB = -L${FFTW3_ROOT}/lib64 -lfftw3

ALL_EXTERNAL_LIBS = $(NETCDF_LIB) $(FFTW_LIB) $(BLAS_LIB) $(LAPACK_LIB) $(XTRA_LIBS)
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
ifeq ($(VENDOR), gcc)
	-@echo  "mpif90 details:       " `mpif90 -v >& tmp; cat tmp | grep GCC; rm tmp`
endif
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
	-@echo  " "	
	-@echo  "--------------------------------------------------------------------------------------"
