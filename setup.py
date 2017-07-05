from setuptools import setup

setup (
  name="fstd2nc",
  version="0.20170705",
  description = 'Converts FSTD files (from Environment Canada) to netCDF files.',
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
  py_modules = ['fstd2nc'],
  install_requires = ['numpy','pytz','netcdf4'],
  package_data = {}, #TODO: language files
  entry_points={
    'console_scripts': [
      'fstd2nc = fstd2nc:_fstd2nc_cmdline_trapped',
    ],
  },

)

