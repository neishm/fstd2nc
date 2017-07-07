# Shortcut for running fstd2nc_locale Makefile

all:
	$(MAKE) -C fstd2nc_locale
%:
	$(MAKE) -C fstd2nc_locale $@

