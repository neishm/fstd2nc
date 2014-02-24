#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

# Extension modules
fstd_core = Extension ('pygeode.formats.fstd_core', sources=['fstd_core.c'], libraries=['rmn'])
fstd_extern = Extension ('pygeode.formats.fstd_extern', sources=['fstd_externmodule.c','fortranobject.c'], extra_objects=['set_igs.o'], libraries=['rmn'])

# PyGeode installation script

setup (	name="python-pygeode-rpn",
	version="0.8.3",
        author="Mike Neish",
	ext_modules=[fstd_core, fstd_extern],
        include_dirs = [np.get_include()],
	packages=["pygeode.plugins.rpn","pygeode.formats"]
)

