Overview
========
This module provides a mechanism for converting between FSTD and netCDF file formats, either through Python or the command-line.

Basic Usage
===========

From the command-line
---------------------
```
python -m fstd2nc [options] <infile> <outfile>
```

From Python
-----------
```python
import fstd2nc
data = fstd2nc.Buffer()
data.read_fstd_file("myfile.fst")
data.write_nc_file("myfile.nc")
```

Requirements
============
This package requires [Python-RPN](https://github.com/meteokid/python-rpn) for reading/writing FSTD files, and [netcdf4-python](https://github.com/Unidata/netcdf4-python) for reading/writing netCDF files.

