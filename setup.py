#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

# Extension modules
fstd_core = Extension ('pygeode.formats.fstd_core', sources=['fstd_core.c'], libraries=['rmn'])
fstd_extern = Extension ('pygeode.formats.fstd_extern', sources=['set_igs.F90'], libraries=['rmn'])

# PyGeode installation script

setup (	name="python-pygeode-rpn",
	version="0.7.1",
        author="Mike Neish",
	ext_modules=[fstd_core, fstd_extern],
	packages=["pygeode.plugins.rpn"]
)

