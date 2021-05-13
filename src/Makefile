#=====================================================	
#  FlowSolve Makefile for src directory
#=====================================================	
SHELL=/bin/sh
include Make.inc

LIBS= ${ALL_EXTERNAL_LIBS}

FS_SRC = data_modules.f90 user_params_module.f90 fourier_differentiation_tools.f90 \
		 BernoulliCosine.f90 solver_utilities.f90 transpose_routines.f90   \
         preliminary_tasks.f90 allocate_memory.f90 mpi_utils.f90 decomposition2D.f90  \
         read_userdata.f90 flow_solve_user_routines.f90   \
         initial_conditions.f90 write_results.f90 user_analysis.f90 apply_forcing.f90 \
         explicit_rhs.f90 explicit_time_step.f90 boundary_routines.f90 transform_xy.f90  \
         pressure_projection.f90 diffusion.f90 flow_solve.f90

FS_OBJS = data_modules.o user_params_module.o fourier_differentiation_tools.o    \
		BernoulliCosine.o solver_utilities.o transpose_routines.o preliminary_tasks.o     \
		allocate_memory.o mpi_utils.o decomposition2D.o read_userdata.o   \
		flow_solve_user_routines.o initial_conditions.o  write_results.o  \
		user_analysis.o explicit_rhs.o apply_forcing.o explicit_time_step.o boundary_routines.o  \
		transform_xy.o pressure_projection.o diffusion.o  module.o setup.o comm2D.o flow_solve.o
		
DRIVER_SRCS = test_driver.f90 data_modules.f90 
			
DRIVER_OBJS = test_driver.o data_modules.o 


#==================
# Targets
#==================

all:
	@echo
	make flow.x
	mv   flow.x ${BASE}	
	
flow.x:	$(FS_SRCS) $(FS_OBJS)
	${FLINKER} ${LDFLAGS} -o $@ ${FS_OBJS} $(LIBS)

test_driver.x: ${DRIVER_OBJS}
	 ${FLINKER} ${LDFLAGS} -o $@ ${DRIVER_OBJS} $(LIBS)
	        
clean:
	@echo
	@echo "src:  removing  *.o, *.mod,  Make.inc, flow_solve_user.f90 "
	rm -f *.o *.mod  Make.inc flow_solve_user.f90 *.x fort*

#=============================================================
# Dependencies  
#=============================================================
flow.x:	$(FS_OBJS)





#=============================================================
# Suffix Rules
# Gfortran creates a file name.mod for each module name. 
# This confuses the gnu make, which thinks that .mod files
# are Modula2 source code files. Override w/ %.o : %.mod
#=============================================================

.SUFFIXES: $(SUFFIXES) .f90  .F90 .f .mod

.f.o:
	${F77} $(DFLAGS) $(F77OPTFLAGS) -c $<

.f90.o:
	${F90} $(DFLAGS) $(FOPTFLAGS) $(ALL_INCLUDES) -c $<
	

.f90.mod:
	${F90} $(FOPTFLAGS) $(ALL_INCLUDES) -c $<
	

%.o : %.mod
