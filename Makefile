include Make.inc
#--------------------------------------------------------------------
#  root directory for "output"
#    i.e. you can specify a remote directory where the output data
#    will actually be written, a symbolic link will be created to
#    a local "output" directory
#--------------------------------------------------------------------
OUTPUT_ROOT = ${BASE}
# either way, don't change this next line
OUTPUT = ${OUTPUT_ROOT}/output
	
all:  
	@echo
	make clean
	make show	
	make test_driver.x 
	make flow.x
	make outdirs
	@echo "Done."

outdirs:
	@echo
	@echo "REMOVING EXISTING OUTPUT DIRECTORIES"
	@echo
	if [ -e "$(OUTPUT_ROOT)" ]; then \
		cd ${OUTPUT_ROOT} ; \
		rm -rf output; \
	fi
	@echo
	@echo "Deleting any local link to output directories"
	@echo
	-rm -rf output
	@echo
	@echo "Creating output root directory if necessary"
	@echo
	if [ ! -e "$(OUTPUT_ROOT)" ]; then \
		mkdir $(OUTPUT_ROOT); \
	fi
	@echo
	@echo "... creating new output directory tree: $(OUTPUT)"
	@echo
	cd ${OUTPUT_ROOT} ; mkdir output output/1D output/2D output/3D output/debug_data
	cd ${OUTPUT_ROOT} ; mkdir output/TS output/codes_etc output/codes_etc/input
	cd ${OUTPUT_ROOT} ; mkdir output/slices output/slices/1D output/slices/2D output/slices/3D output/movie_frames
	cd ${BASE} ; cp Make.inc ${OUTPUT_ROOT}/output/codes_etc/.
	if [ -e "${BASE}/FlowSolveConfig" ]; then \
		cp ${BASE}/FlowSolveConfig ${OUTPUT_ROOT}/output/codes_etc/. ; \
	fi
	cd ${BASE}/input ; cp -r * ${OUTPUT_ROOT}/output/codes_etc/input/.
	cd ${BASE} ; cp -r src ${OUTPUT_ROOT}/output/codes_etc/.
	cd ${OUTPUT_ROOT}/output/codes_etc/src ; rm -f *.o *.mod *.x
	cd ${OUTPUT_ROOT}/output/codes_etc/src/p3dfft ; rm -f *.o *.mod *.a
	cd ${OUTPUT_ROOT} ; rm -rf output/codes_etc/input/*pdf output/codes_etc/input/*eps output/codes_etc/input/*png output/codes_etc/input/*nc
	@echo
	@echo "... if needed, creating local link to output directory: $(OUTPUT)"
	@echo
	@echo "... BASE $(BASE)"	
	if [ ! -e "${BASE}/output" ]; then \
		cd ${BASE} ; \
		ln -s ${OUTPUT} output ; \
	fi	
	@echo
	cd ${BASE} ;  make show	
	@echo "Done."

				
flow.x:
	@echo
	@echo "... building executable flow.x"
	@echo "... "
	@echo "... compiling 2d transpose routines"
	cp -f Make.inc src/p3dfft
	cd src/p3dfft; make "BASE=${BASE}" dot_o_files
	@echo "... "
	@echo "... making temp directory for compilation"
	rm -rf compile_directory 
	mkdir compile_directory
	@echo "... "
	@echo "... copying original source code to temp directory"
	cp -rf Make.inc Makefile src compile_directory/.
	cp -f src/p3dfft/*.o src/p3dfft/*.mod compile_directory/src/. 
	@echo "... "
	@echo "... copying standard user routines to use for compilation"
	-cp -f input/*.f90  input/Makefile compile_directory/src/. # >& /dev/null
	@echo "... "
	@echo "... compiling flow.x in compile_directory/src directory"
	cp Make.inc compile_directory/src/. ; cd compile_directory/src; make "BASE=${BASE}" all
	@echo "... "
	@echo "... removing temp directory. cleaning up"
	rm -rf compile_directory   
	rm -f FlowSolveConfig
	make showconfig > FlowSolveConfig
	@echo "Done."
	@echo
		
		
test_driver.x:	src/test_driver.f90 
	@echo
	@echo "... building executable test_driver.x"
	rm -f ${BASE}/test_driver*; clear
	rm -f output/*vec; clear
	cp -f Make.inc src; cd src; make showconfig; make "BASE=${BASE}" test_driver.x
	mv src/test_driver.x ${BASE}/test_driver.x
	@echo "test_driver.x built. It can be run using mpirun:  e.g.  mpirun -np 5 ./test_driver.x"
	@echo
		
		
clean:
	@echo
	make cleansrcdir
	rm -rf compile_directory
	rm -f flow.x test_driver* install_test* 
	rm -f runlog FlowSolveConfig 
	rm -f svn_status fort.* *.tmp a.out *.aux *.log *.synctex.gz *.pdf *.toc
	rm -f *~
	rm -f input/data_tools/python_scripts/*.pyc input/data_tools/*.pyc input/data_tools/paths/*
	clear
	@echo " "
	@echo "---------------- "
	@echo "Top Directory"
	@echo "---------------- "
	ls 
	@echo " "
	@echo "---------------- "
	@echo "Source Directory"
	@echo "---------------- "
	ls src
	@echo " "
	@echo "---------------- "
	@echo "Input Directory"
	@echo "---------------- "
	if [ -e "input" ]; then \
	 ls input ; \
	fi
	@echo " "
	@echo "---------------- "
	@echo "Output Directory"
	@echo "---------------- "
	cd output;pwd;ls
	@echo "Done."
	
cleansrcdir:  	
	@echo
	cp -f Make.inc src
	cd src;  make clean
	cp -f Make.inc src/p3dfft
	cd src/p3dfft;  make clean ; rm -f Make.inc
	
show:
	@echo
	rm -f FlowSolveConfig
	make "BASE=${BASE}" -f  Make.inc showconfig >& FlowSolveConfig

help:
	@echo
	@echo "===========  available make targets ==========================="
	@echo "make show             -- show configuration variables"
	@echo "make all              -- clean,outdirs,libs,flow.x"
	@echo "make outdirs          -- create new output directories"
	@echo "make clean            -- cleans all dirs except output"
	@echo "make flow.x           -- create executable flow.x"
	@echo "make test_driver.x    -- create executable to test libraries"
	@echo "Done." 
	@echo "==============================================================="

