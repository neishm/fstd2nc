#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

# Extension modules
fstd_core = Extension ('pygeode.formats.fstd_core', sources=['fstd_core.c'], libraries=['rmn'])

# PyGeode installation script

setup (	name="python-pygeode-rpn",
	version="2.x-UNRELEASED",
        author="Mike Neish",
	ext_modules=[fstd_core],
        include_dirs = [np.get_include()],
	packages=["pygeode.formats"],
        package_dir = {"pygeode.formats":""},
	scripts=["fstd2nc"]
)

