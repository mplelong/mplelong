include Make.inc
FF = ${F90}
FFLAGS = ${F77OPTFLAGS}

INCL= 
LDFLAGS=
DFLAGS = "-DDOUBLE_PREC"


P3DFFT_ROOT = .
P3DFFT_LIB = ${P3DFFT_ROOT}/libp3dfft.3.a


FFT3DLIB = module.o setup.o comm2D.o

 
dot_o_files: $(FFT3DLIB)
	
test: $(FFT3DLIB) test1_2D.o
	$(FF) test1_2D.o -L. -lp3dfft.3 $(LDFLAGS) 

setup.o: setup.f module.o
comm2D.o: comm2D.F module.o

.c.o: 
		$(CC) -c $(CFLAGS) *.c
.F.o:
		$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o: 
		$(FF) -c $(DFLAGS) $(FFLAGS) $(INCL) $<
clean:
	/bin/rm -f $(FFT3DLIB) *.mod *.i a.out test1_2D.o *~
	
	
