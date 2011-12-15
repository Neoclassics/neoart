include ../src/Makefile-neoart.include

TEST_OBJECTS = test.o c1.o c2.o c3.o k22.o psflux.o \
	TEST1.o TEST2.o TEST3.o TEST4.o TEST5.o\
	TEST6.o TEST7.o TEST8.o TEST9.o TEST10.o TEST11.o \
	TEST12.o TEST13.o TEST14.o TEST16.o TEST18.o TEST19.o

all: test run_test

test: $(addprefix ../src/, $(NEOART_OBJECTS)) $(TEST_OBJECTS)
	$(FC) $(FFLAGS) -o test $^

run_test:
	./test > test_output
	diff test_output test_output_reference

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f90 
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o test


