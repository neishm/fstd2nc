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
  --progress            Display a progress bar during the conversion. Requires
                        the "progress" module.
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

Using in a Python script
========================

Simple conversion
--------------------------------------
```python
import fstd2nc
data = fstd2nc.Buffer("myfile.fst")
data.write_nc_file("myfile.nc")
```

You can control `fstd2nc.Buffer` using parameters similar to the command-line arguments.  The usual convention is *--arg-name* from the command-line would be passed as *arg_name* from Python.

For example:
```python
import fstd2nc
# Strip out FSTD metadata from the dataset, and select only TT,HU variables.
data = fstd2nc.Buffer("myfile.fst", minimal_metadata=True, vars=['TT','HU'])
# Set the reference date to Jan 1, 2000 in the netCDF file.
data.write_nc_file("myfile.nc", reference_date='2000-01-01')
```

Interfacing with xarray
---------------------------------------------------------------------------------

For more complicated conversions, you can manipulate the data as an `xarray.Dataset` object:
```python
import fstd2nc

# Open the FSTD file, and use dates of validity for the time axis.
data = fstd2nc.Buffer("myfile.fst", squash_forecasts=True)

# Access the data as an xarray.Dataset object.
dataset = data.to_xarray()
print (dataset)

# Convert surface pressure to Pa.
dataset['P0'] *= 100
dataset['P0'].attrs['units'] = 'Pa'

# (Can further manipulate the dataset here)
# ...

# Write the final result to netCDF using xarray:
dataset.to_netcdf("myfile.nc")
```

Installing
==========

The easiest way to install is using [pip](https://pip.pypa.io/en/stable):
```
pip install fstd2nc
```

If you're processing many input files into a single netCDF file, you could get some useful features (progress bar, quick file scans) by running:
```
pip install fstd2nc[manyfiles]
```


Using in a Pydap server
=======================

This package includes a handler for [Pydap](https://github.com/pydap/pydap), which enables you to serve your FSTD files via the OPeNDAP protocol.

To install all the pre-requisites:
```
pip install pydap fstd2nc[dap]
```

You can then test it by running
```
pydap -d [your data directory]
```


Requirements
============

Basic requirements
--------------------

This package requires [Python-RPN](https://github.com/meteokid/python-rpn) for reading/writing FSTD files, and [netcdf4-python](https://github.com/Unidata/netcdf4-python) for reading/writing netCDF files.

Optional requirements
---------------------

For reading large numbers of input files (>100), this utility can leverage [pandas](https://github.com/pandas-dev/pandas) to quickly process the FSTD record headers.

The [progress](https://github.com/verigak/progress) module is required in order to use the `--progress` option.

The `.to_xarray()` Python method requires the [xarray](https://github.com/pydata/xarray) and [dask](https://github.com/dask/dask) packages.

