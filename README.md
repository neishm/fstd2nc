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
  --no-progress         Disable the progress bar.
  --minimal-metadata    Don't include RPN record attributes and other internal
                        information in the output metadata. This is the
                        default behaviour.
  --rpnstd-metadata     Include all RPN record attributes in the output
                        metadata.
  --rpnstd-metadata-list nomvar,...
                        Specify a minimal set of RPN record attributes to
                        include in the output file.
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
  --datev, --squash-forecasts
                        Use the date of validity for the "time" axis. This is
                        the default.
  --dateo, --forecast-axis
                        Use the date of original analysis for the time axis,
                        and put the forecast times into a separate "forecast"
                        axis.
  --ensembles           Collect different etikets for the same variable
                        together into an "ensemble" axis.
  --profile-momentum-vars VAR1,VAR2,...
                        Comma-separated list of variables that use momentum
                        levels.
  --profile-thermodynamic-vars VAR1,VAR2,...
                        Comma-separated list of variables that use
                        thermodynamic levels.
  --missing-bottom-profile-level
                        Assume the bottom level of the profile data is
                        missing.
  --strict-vcoord-match
                        Require the IP1/IP2/IP3 parameters of the vertical
                        coordinate to match the IG1/IG2/IG3 paramters of the
                        field in order to be used. The default behaviour is to
                        use the vertical record anyway if it's the only one in
                        the file.
  --diag-as-model-level
                        Treat diagnostic (near-surface) data as model level
                        '1.0'. This is the default behaviour.
  --split-diag-level    Put the diagnostic (near-surface) data in a separate
                        variable, away from the 3D model output. Suffices will
                        be added to distinguish the different types of levels
                        (e.g. _vgrid4 and _vgrid5 for diagnostic height and
                        hybrid levels respectively).
  --ignore-diag-level   Ignore data on diagnostic (near-surface) height.
  --subgrid-axis        For data on supergrids, split the subgrids along a
                        "subgrid" axis. The default is to leave the subgrids
                        stacked together as they are in the RPN file.
  --keep-LA-LO          Include LA and LO records, even if they appear to be
                        redundant.
  --filter CONDITION    Subset RPN file records using the given criteria. For
                        example, to convert only 24-hour forecasts you could
                        use --filter ip2==24
  --exclude NAME,NAME,...
                        Exclude some axes or derived variables from the
                        output. Note that axes will only be excluded if they
                        have a length of 1.
  --metadata-file METADATA_FILE
                        Use metadata from the specified file. You can repeat
                        this option multiple times to build metadata from
                        different sources.
  --rename OLDNAME=NEWNAME,...
                        Apply the specified name changes to the variables.
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
  --compression COMPRESSION
                        Compression level for the netCDF file. Only used if
                        --zlib is set. Default: 4.
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
data.to_netcdf("myfile.nc")
```

You can control `fstd2nc.Buffer` using parameters similar to the command-line arguments.  The usual convention is *--arg-name* from the command-line would be passed as *arg_name* from Python.

For example:
```python
import fstd2nc
# Select only TT,HU variables.
data = fstd2nc.Buffer("myfile.fst", vars=['TT','HU'])
# Set the reference date to Jan 1, 2000 in the netCDF file.
data.to_netcdf("myfile.nc", reference_date='2000-01-01')
```

Interfacing with xarray
---------------------------------------------------------------------------------

For more complicated conversions, you can manipulate the data as an [xarray.Dataset](http://xarray.pydata.org/en/stable/data-structures.html#dataset) object by using the `to_xarray()` method:
```python
import fstd2nc

# Open the FSTD file.
data = fstd2nc.Buffer("myfile.fst")

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

Interfacing with iris
---------------------------------------------------------------------------------

You can interface with [iris](https://scitools.org.uk/iris/docs/latest/index.html) by using the `.to_iris()` method (requires iris version 2.0 or greater).
This will give you an [iris.cube.CubeList](https://scitools.org.uk/iris/docs/latest/iris/iris/cube.html#iris.cube.CubeList) object:
```python
import fstd2nc
import iris.quickplot as qp
from matplotlib import pyplot as pl

# Open the FSTD file.
data = fstd2nc.Buffer("myfile.fst")

# Access the data as an iris.cube.CubeList object.
cubes = data.to_iris()
print (cubes)

# Plot all the data (assuming we have 2D fields)
for cube in cubes:
  qp.contourf(cube)
  pl.gca().coastlines()

pl.show()
```

Interfacing with pygeode
---------------------------------------------------------------------------------

You can create a [pygeode.Dataset](http://pygeode.github.io/dataset.html) object using the `.to_pygeode()` method (requires pygeode version 1.2.2 or greater):
```python
import fstd2nc

# Open the FSTD file.
data = fstd2nc.Buffer("myfile.fst")

# Access the data as a pygeode.Dataset object.
dataset = data.to_pygeode()
print (dataset)
```


Installing
==========

The easiest way to install is using [pip](https://pip.pypa.io/en/stable):
```
pip install fstd2nc
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

The `.to_xarray()` Python method requires the [xarray](https://github.com/pydata/xarray) and [dask](https://github.com/dask/dask) packages.

The `.to_iris()` Python method requires the [iris](https://scitools.org.uk/iris/docs/latest/index.html) package, along with the `.to_xarray()` dependencies.

The `.to_pygeode()` Python method requires the [pygeode](https://github.com/pygeode/pygeode) package, along with the `.to_xarray()` dependencies.
