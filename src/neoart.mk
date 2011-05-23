

FC = mpif90 
SRCDIR = ../src
RUNDIR = ../run
TESTDIR = ../tst
FFLAGS = -fdefault-real-8 

# -O3 -ftracer -fomit-frame-pointer -pipe -fweb -fdefault-real-8


vpath %.f $(SRCDIR)

exec :  RUNME.o c1.o c2.o c3.o k22.o psflux.o ralf.o    \
	dandv.o ps.o bp.o viscos.o visfus.o viscol.o get_geom.o    \
	menn.o penq.o colxi.o perr.o geom.o circgeom.o visgeom.o   \
	erf.o ludcmp.o lubksb.o advance.o interp.o neoart.o \
	class.o
	$(FC) $(FFLAGS) -o $(RUNDIR)/neoart RUNME.o c1.o c2.o c3.o k22.o psflux.o ralf.o   \
	dandv.o ps.o bp.o viscos.o visfus.o viscol.o  get_geom.o neoart.o   \
	menn.o penq.o colxi.o perr.o geom.o circgeom.o visgeom.o \
	erf.o ludcmp.o lubksb.o advance.o interp.o class.o 


%.o : %.f
	$(FC) $(FFLAGS) -c -o $@ $<

test :  test.o c1.o c2.o c3.o k22.o psflux.o ralf.o    \
	dandv.o ps.o bp.o viscos.o visfus.o viscol.o get_geom.o    \
	menn.o penq.o colxi.o perr.o geom.o circgeom.o visgeom.o   \
	erf.o ludcmp.o lubksb.o advance.o interp.o neoart.o \
	class.o TEST1.o 
	$(FC) $(FFLAGS) -o $(TESTDIR)/test test.o c1.o c2.o c3.o k22.o psflux.o ralf.o   \
	dandv.o ps.o bp.o viscos.o visfus.o viscol.o  get_geom.o neoart.o   \
	menn.o penq.o colxi.o perr.o geom.o circgeom.o visgeom.o \
	erf.o ludcmp.o lubksb.o advance.o interp.o class.o TEST1.o

