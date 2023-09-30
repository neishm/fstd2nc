###############################################################################
# Copyright 2017-2023 - Climate Research Division
#                       Environment and Climate Change Canada
#
# This file is part of the "fstd2nc" package.
#
# "fstd2nc" is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# "fstd2nc" is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with "fstd2nc".  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

from setuptools import setup, find_packages
from fstd2nc import __version__

with open("README.md","r") as f:
  long_description = f.read()

setup (
  name="fstd2nc",
  version=__version__,
  description = 'Converts RPN standard files (from Environment Canada) to netCDF files.',
  long_description = long_description,
  # https://stackoverflow.com/a/26737258/9947646
  long_description_content_type='text/markdown',
  url = 'https://github.com/neishm/fstd2nc',
  author="Mike Neish",
  license = 'LGPL-3',
  classifiers = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
  ],
  packages = find_packages(),
  setup_requires = ['pip >= 8.1'],
  install_requires = ['numpy >= 1.13.0, != 1.15.3','netcdf4','fstd2nc-deps >= 0.20200304.0','progress'],
  extras_require = {
    'manyfiles': ['pandas'],
    'array': ['xarray>=0.10.3','dask','toolz'],
    'iris': ['iris>=2.0','xarray>=0.10.3','dask','toolz'],
    'pygeode': ['pygeode>=1.2.2','xarray>=0.10.3','dask','toolz'],
  },
  package_data = {
    'fstd2nc': ['locale/*/LC_MESSAGES/fstd2nc.mo'],
  },
  entry_points={
    'console_scripts': [
      'fstd2nc = fstd2nc.__main__:_fstd2nc_cmdline_trapped',
      'fstdump = fstd2nc.__main__:_fstdump',
      'ccc2nc = cccbuffer.__main__:run',
      'cccdump = cccbuffer.cccdump:run',
    ],
    'xarray.backends': [
      'fstd = fstd2nc._xarray_hook:FSTDBackendEntrypoint',
      'ccc = cccbuffer._xarray_hook:CCCBackendEntrypoint',
    ],
  },

)

