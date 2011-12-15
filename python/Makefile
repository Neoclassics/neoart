include ../src/Makefile-neoart.include

neoart_so = neoart/_neoart.so

all: $(neoart_so)

tests: $(neoart_so)
	nosetests -w tests

$(neoart_so):
	python setup.py config_fc --f77flags="$(FFLAGS)"\
		--f90flags="$(FFLAGS)"\
		build_ext --inplace

clean:
	rm -rf build src/_neoart-f2pywrappers.f src/_neoartmodule.c
	rm -f $(neoart_so) neoart/*.pyc

.PHONY: tests
