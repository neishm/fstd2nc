###############################################################################
# Copyright 2017 - Climate Research Division
#                  Environment and Climate Change Canada
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

setup (
  name="fstd2nc",
  version=__version__,
  description = 'Converts RPN standard files (from Environment Canada) to netCDF files.',
  url = 'https://github.com/neishm/fstd2nc',
  author="Mike Neish",
  license = 'LGPL-3',
  classifiers = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
  ],
  packages = find_packages(),
  py_modules = ['fstd2nc','fstd2dap','fstd2nc_extra'],
  setup_requires = ['pip >= 8.1'],
  install_requires = ['numpy >= 1.13.0','netcdf4','fstd2nc-deps'],
  extras_require = {
    'dap': ['Pydap[server,functions]'],
  },
  package_data = {
    'fstd2nc_locale': ['*/LC_MESSAGES/fstd2nc.mo'],
  },
  entry_points={
    'console_scripts': [
      'fstd2nc = fstd2nc:_fstd2nc_cmdline_trapped',
    ],
    'pydap.handler': [
      'fstd = fstd2dap:FST_Handler',
    ],
  },

)

