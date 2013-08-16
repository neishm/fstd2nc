SUBDIRS = 
LIBS = rpn bmf geophy
COMMON = io compress
include rules.mk

install:
	$(MAKE)
	gfortran -fPIC -c set_igs.F90
	f2py set_igs.F90 -m fstd_extern
	cp /usr/share/pyshared/numpy/f2py/src/fortranobject.[ch] .
	python setup.py build
	python setup.py install --prefix=$(DESTDIR)/usr/local
	./cp_libs.sh
	./cp_bins.sh
