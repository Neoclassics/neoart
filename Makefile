all : neoart test

neoart :
	gmake -C obj -f ../src/neoart.mk exec

test: 
	gmake -C obj -f ../src/neoart.mk test
	tst/test > tst/test_output 
	diff tst/test_output tst/test_output_reference

clean : 
	rm obj/*.o
