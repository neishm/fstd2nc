# Shortcut for running fstd2nc_locale Makefile

all:
	$(MAKE) -C fstd2nc/locale
%:
	$(MAKE) -C fstd2nc/locale $@

