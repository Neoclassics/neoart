include Makefile-neoart.include

neoart: $(NEOART_OBJECTS) RUNME.o
	$(FC) $(FFLAGS) -o ../neoart $^

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f90 
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o ../neoart

