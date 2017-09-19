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
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --minimal-metadata    Don't include RPN record attributes and other internal
                        information in the output metadata.
  --ignore-typvar       Tells the converter to ignore the typvar when deciding
                        if two records are part of the same field. Default is
                        to split the variable on different typvars.
  --ignore-etiket       Tells the converter to ignore the etiket when deciding
                        if two records are part of the same field. Default is
                        to split the variable on different etikets.
  --vars VAR1,VAR2,...  Comma-separated list of variables to convert. By
                        default, all variables are converted.
  --fill-value FILL_VALUE
                        The fill value to use for masked (missing) data. Gets
                        stored as '_FillValue' attribute in the metadata.
                        Default is '1e+30'.
  --squash-forecasts    Use the date of validity for the "time" axis.
                        Otherwise, the default is to use the date of original
                        analysis, and the forecast length goes in a "forecast"
                        axis.
  --subgrid-axis        For data on supergrids, split the subgrids along a
                        "subgrid" axis. The default is to leave the subgrids
                        stacked together as they are in the RPN file.
  --filter CONDITION    Subset RPN file records using the given criteria. For
                        example, to convert only 24-hour forecasts you could
                        use --filter ip2==24
  --metadata-file METADATA_FILE
                        Use metadata from the specified file. You can repeat
                        this option multiple times to build metadata from
                        different sources.
  --time-units {seconds,minutes,hours,days}
                        The units for the output time axis. Default is hours.
  --reference-date YYYY-MM-DD
                        The reference date for the output time axis. The
                        default is the starting date in the RPN file.
  --msglvl {0,DEBUG,2,INFORM,4,WARNIN,6,ERRORS,8,FATALE,10,SYSTEM,CATAST}
                        How much information to print to stdout during the
                        conversion. Default is WARNIN.
  --nc-format {NETCDF4,NETCDF4_CLASSIC,NETCDF3_CLASSIC,NETCDF3_64BIT_OFFSET,NETCDF3_64BIT_DATA}
                        Which variant of netCDF to write. Default is NETCDF4.
  --zlib                Turn on compression for the netCDF file. Only works
                        for NETCDF4 and NETCDF4_CLASSIC formats.
  -f, --force           Overwrite the output file if it already exists.
  --no-history          Don't put the command-line invocation in the netCDF
                        metadata.
```

From Python
-----------
Simple conversion to netCDF:
```python
import fstd2nc
data = fstd2nc.Buffer("myfile.fst")
data.write_nc_file("myfile.nc")
```
To use the higher-level data structures in Python for other purposes:
```python
import fstd2nc
import numpy as np
data = fstd2nc.Buffer("myfile.fst")

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

