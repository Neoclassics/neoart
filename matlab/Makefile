include ../src/Makefile-neoart.include

FOPTIMFLAGS=-g

all: neoart.mex*

neoart.mex*: $(addprefix ../src/, $(NEOART_OBJECTS)) neoart_top.f90
	mex -fortran 'FC=${FC}' 'FOPTIMFLAGS=${FOPTIMFLAGS}' -output neoart $^

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f90 
	$(FC) $(FFLAGS) -c -o $@ $<

clean: 
	rm -f neoart.mex*
