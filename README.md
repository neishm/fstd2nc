[Version française](README_fr.md)

Overview
========
This module provides a mechanism for converting between FSTD and netCDF file formats, either through Python or the command-line.


Installing
==========

The easiest way to install is using [pip](https://pip.pypa.io/en/stable):
```
pip install fstd2nc
```


Basic Usage
===========

From the command-line
---------------------
```
python -m fstd2nc [options] <infile(s)> <outfile>

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --no-progress         Disable the progress bar.
  --serial              Disables multithreading/multiprocessing. Useful for
                        resource-limited machines.
  --minimal-metadata    Don't include internal record attributes and other
                        internal information in the output metadata. This is
                        the default behaviour.
  --internal-metadata, --rpnstd-metadata
                        Include all internal record attributes in the output
                        metadata.
  --metadata-list nomvar,..., --rpnstd-metadata-list nomvar,...
                        Specify a minimal set of internal record attributes to
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
  --accum-vars NAME,NAME,...
                        Specify which fields to treat as accumulated
                        quantities (using IP3 as accumulation period).
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
  --vardict VARDICT     Use metadata from the specified variable dictionary
                        (XML format).
  --opdict              Similar to above, but use the standard CMC-RPN
                        operational dictionary.
  --sfc-agg-vars NAME,NAME,...
                        Define additional surface aggregate fields.
  --soil-depths SOIL_DEPTHS
                        Define custom depths for soil fields (WSOL,ISOL).
                        Defaults are 0.05,0.1,0.2,0.4,1.0,2.0,3.0.
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
                        (i.e. _diag_level and _model_levels for diagnostic
                        height and hybrid levels respectively).
  --ignore-diag-level   Ignore data on diagnostic (near-surface) height.
  --only-diag-level     Only use the diagnostic (near-surface) height,
                        ignoring other atmospheric levels.
  --thermodynamic-levels, --tlev
                        Only convert data that's on 'thermodynamic' vertical
                        levels.
  --momentum-levels, --mlev
                        Only convert data that's on 'momentum' vertical
                        levels.
  --vertical-velocity-levels, --wlev
                        Only convert data that's on 'vertical velocity'
                        levels.
  --subgrid-axis        For data on supergrids, split the subgrids along a
                        "subgrid" axis. The default is to leave the subgrids
                        stacked together as they are in the RPN file.
  --keep-LA-LO          Include LA and LO records, even if they appear to be
                        redundant.
  --no-adjust-rlon      For rotated grids, do NOT adjust rlon coordinate to
                        keep the range in -180..180. Allow the rlon value to
                        be whatever librmn says it should be.
  --bounds              Include grid cell boundaries in the output.
  --filter CONDITION    Subset RPN standard file records using the given
                        criteria. For example, to convert only 24-hour
                        forecasts you could use --filter ip2==24
  --exclude NAME,NAME,...
                        Exclude some axes, attributes, or derived variables
                        from the output. For instance, excluding
                        'leadtime,reftime' can help for netCDF tools that
                        don't recognize leadtime and reftime as valid
                        coordinates. Note that axes will only be excluded if
                        they have a length of 1.
  --yin                 Select first subgrid from a supergrid.
  --yang                Select second subgrid from a supergrid.
  --crop-to-smallest-grid
                        Crop grids to the smaller (inner core) domain for LAM
                        outputs.
  --metadata-file METADATA_FILE
                        Use metadata from the specified file. You can repeat
                        this option multiple times to build metadata from
                        different sources.
  --rename OLDNAME=NEWNAME,...
                        Apply the specified name changes to the variables.
  --conventions CONVENTIONS
                        Set the "Conventions" attribute for the netCDF file.
                        Default is "CF-1.6". Note that this has no effect on
                        the structure of the file.
  --no-conventions      Omit the "Conventions" attribute from the netCDF file
                        entirely. This can help for netCDF tools that have
                        trouble recognizing the CF conventions encoded in the
                        file.
  --time-units {seconds,minutes,hours,days}
                        The units for the output time axis. Default is hours.
  --reference-date YYYY-MM-DD
                        The reference date for the output time axis. The
                        default is the starting date in the RPN standard file.
  --fstd-compat         Adds a compatibility layer to the netCDF output file,
                        so it can also function as a valid FSTD file.
                        EXPERIMENTAL.
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
  -q, --quiet           Don't display any information except for critical
                        error messages. Implies --no-progress.
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

You can control `fstd2nc.Buffer` using parameters similar to the command-line arguments.  The convention is that *--arg-name* from the command-line would be passed as *arg_name* from Python.

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

For more complicated conversions, you can manipulate the data as an [xarray.Dataset](http://xarray.pydata.org/en/stable/data-structures.html#dataset) object:
```python
import xarray as xr

# Open the FSTD file.
# Access the data as an xarray.Dataset object.
dataset = xr.open_dataset("myfile.fst", engine="fstd")
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

Interfacing with fstpy
---------------------------------------------------------------------------------

You can load data from an [fstpy](https://gitlab.science.gc.ca/CMDS/fstpy) table using the `.from_fstpy()` method (requires fstpy version 2.1.9 or greater):
```python
import fstd2nc
import fstpy
table = fstpy.StandardFileReader('myfile.fst').to_pandas()
data = fstd2nc.Buffer.from_fstpy(table)
```
You can also export to an fstpy table using the `.to_fstpy()` method:
```python
import fstd2nc
table = fstd2nc.Buffer('myfile.fst').to_fstpy()
```


Requirements
============

Basic requirements
--------------------

This package requires [Python-RPN](https://github.com/meteokid/python-rpn) for reading/writing FSTD files, and [netcdf4-python](https://github.com/Unidata/netcdf4-python) for reading/writing netCDF files.

Optional dependencies
---------------------

A useful variable dictionary for the `--vardict` option is available [here](https://collaboration.cmc.ec.gc.ca/cmc/CMOI/VariableDictionary/).

For reading large numbers of input files (>100), this utility can leverage [pandas](https://github.com/pandas-dev/pandas) to quickly process the FSTD record headers.

The `.to_xarray()` Python method requires the [xarray](https://github.com/pydata/xarray) and [dask](https://github.com/dask/dask) packages.

The `.to_iris()` Python method requires the [iris](https://scitools.org.uk/iris/docs/latest/index.html) package, along with the `.to_xarray()` dependencies.

The `.to_pygeode()` Python method requires the [pygeode](https://github.com/pygeode/pygeode) package, along with the `.to_xarray()` dependencies.

The `.to_fstpy()` Python method requires the [fstpy](https://gitlab.science.gc.ca/CMDS/fstpy) package (*internal link*).
