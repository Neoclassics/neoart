all : run/neoart
	gmake -C obj -f ../src/neoart.mk

clean : 
	rm obj/*.o
