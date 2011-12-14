FC = gfortran
FFLAGS = -fdefault-real-8 -Wuninitialized
FFLAGS += -fPIC # needed to compile matlab interface

export FC
export FFLAGS

all : neoart test

neoart:
	$(MAKE) -C src neoart

matlab:
	$(MAKE) -C matlab

test: tests
tests:
	$(MAKE) -C tests
	(cd tests && ./test > test_output)
	(cd tests && diff test_output test_output_reference)

clean:
	@for y in src tests matlab; do $(MAKE) -C $$y/ clean ; done

.PHONY: matlab test tests

