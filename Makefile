build:
	python setup.py build

install:
	python setup.py install --prefix=$(DESTDIR)/usr/local
	./cp_bins.sh

