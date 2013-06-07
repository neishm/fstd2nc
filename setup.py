#!/usr/bin/env python

from distutils.core import setup, Extension

# Extension modules
fstd_core = Extension ('pygeode.formats.fstd_core', sources=['fstd_core.c'], libraries=['rmn'])

# PyGeode installation script

setup (	name="python-pygeode-rpn",
	version="0.7.1",
        author="Mike Neish",
	ext_modules=[fstd_core],
	packages=["pygeode.plugins.rpn"]
)

