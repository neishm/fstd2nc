Overview
========
This module provides a mechanism for converting between FSTD and netCDF file formats, either through Python or the command-line.

Basic Usage
===========

From the command-line
---------------------
```
python -m fstd2nc [options] <infile> <outfile>

optional arguments:
  --minimal-metadata    Don't include FSTD record attributes and other
                        internal information in the output file.
  --ignore-etiket       Tells the converter to ignore the etiket when deciding
                        if two records are part of the same field. Default is
                        to split the variable on different etikets.
  --vars VAR1,VAR2,...  Comma-seperated list of variables to convert. By
                        default, all variables are converted.
  --fill-value FILL_VALUE
                        The fill value to use for masked (missing) data. Gets
                        stored as '_FillValue' attribute in the netCDF file.
                        Default is '1e+30'.
  --squash-forecasts    Use the date of validity for the "time" axis.
                        Otherwise, the default is to use the date of original
                        analysis, and the forecast length goes in a "forecast"
                        axis.
  --metadata-file METADATA_FILE
                        Apply netCDF metadata from the specified file.
  --time-units {seconds,minutes,hours,days}
                        The units of time for the netCDF file. Default is
                        hours.
  --buffer-size BUFFER_SIZE
                        How much data to write at a time (in MBytes). Default
                        is 100.
```

From Python
-----------
Simple conversion to netCDF:
```python
import fstd2nc
data = fstd2nc.Buffer()
data.read_fstd_file("myfile.fst")
data.write_nc_file("myfile.nc")
```
To use the higher-level data structures in Python for other purposes:
```python
import fstd2nc
import numpy as np
data = fstd2nc.Buffer()
data.read_fstd_file("myfile.fst")

for name, atts, axes, array in data:

  # The name of the field
  print (name)

  # 'atts' is a dictionary of metadata associated with this field.
  # Contains info from the FSTD record header(s), as well as extra metadata
  # that's useful for netCDF utilities.
  print (atts)

  # 'axes' is an ordered dictionary containing the coordinate system.
  print (list(axes.keys()))   # Names of the dimensions

  # 'array' points to a "symbolic" array, which isn't yet loaded in memory.
  # That way, you can slice it to select only the subset of data you need.
  # To load the values, pass it to numpy.asarray(), or do a numpy operation
  # on it.
  print (array.shape)
  print (np.mean(array))
```

Requirements
============
This package requires [Python-RPN](https://github.com/meteokid/python-rpn) for reading/writing FSTD files, and [netcdf4-python](https://github.com/Unidata/netcdf4-python) for reading/writing netCDF files.

