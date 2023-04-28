
# Helper function - converts "number of 15-minute intervals since year 1"
# to a datetime object.
from fstd2nc.mixins import vectorize
@vectorize
def year1offset_to_date (timestep):
  from cftime import DatetimeNoLeap
  from datetime import timedelta
  timestep = int(timestep)
  return DatetimeNoLeap(year=1,month=1,day=1) + timedelta(minutes=timestep*15)

# Helper function - converts dates from "yyyymmddhh" format to a datetime
# object.
@vectorize
def yyyymmddhh_to_date (timestep):
  from cftime import DatetimeNoLeap
  from datetime import timedelta
  yyyymmddhh = int(timestep)
  yyyymmdd, hh = divmod(yyyymmddhh,100)
  if hh > 24: raise ValueError
  yyyymm, dd = divmod(yyyymmdd, 100)
  yyyy, mm = divmod(yyyymm, 100)
  return DatetimeNoLeap(year=yyyy,month=mm,day=dd) + timedelta(hours=hh)

# Helper function - converts dates from "yyyymmdd" format to a datetime
# object.
@vectorize
def yyyymmdd_to_date (timestep):
  from cftime import DatetimeNoLeap
  yyyymmdd = int(timestep)
  yyyymm, dd = divmod(yyyymmdd, 100)
  yyyy, mm = divmod(yyyymm, 100)
  return DatetimeNoLeap(year=yyyy,month=mm,day=dd)

# Helper function - converts dates from "yyyymm" format to a datetime
# object.
@vectorize
def yyyymm_to_date (timestep):
  from cftime import DatetimeNoLeap
  yyyymm = int(timestep)
  yyyy, mm = divmod(yyyymm, 100)
  return DatetimeNoLeap(year=yyyy,month=mm,day=1)

# Helper function - extract year / month from filename.
@vectorize
def filename_to_date (filename):
  from cftime import DatetimeNoLeap
  from os.path import basename
  import re
  filename = basename(filename)
  match = re.search("_([0-9]{4})_m([0-9]{2})_", filename)
  if match is None:
    raise ValueError
  yyyy, mm = match.groups()
  return DatetimeNoLeap(year=int(yyyy),month=int(mm),day=1)

# Mixin for decoding dates/times from CCC files.
from fstd2nc.mixins import BufferBase
from fstd2nc.stdout import _, info, warn, error
class Times(BufferBase):
  # This method adds additional command-line arguments needed for the CCC
  # file interpretation.
  @classmethod
  def _cmdline_args (cls, parser):
    super (Times,cls)._cmdline_args(parser)
    parser.add_argument ('--time-type', choices=('year1-offset','yyyymmddhh','yyyymmdd','yyyymm','from-file','code','auto'), help=_('How to interpret the timestep values into a real date/time.'), default='auto')

  # This method is responsible for setting up the initial table of information
  # for the records in the file.
  # It adds any specific metadata columns needed for the mixin.
  # In this case, we are adding/modifying fields related to the CCC file.
  def __init__ (self, *args, **kwargs):
    """
    time_type : str, optional
        How to interpret the timestep values into a real date/time.
        ('year1-offset','yyyymmddhh','yyyymmdd','yyyymm','from-file','code','auto')
    """
    self._time_type = kwargs.pop('time_type','auto')
    super(Times,self).__init__(*args,**kwargs)

  # This routine handles metadata and coordinates for the data
  # (after it has been structured into multi-dimensional arrays).
  # Here we annotate the time axis.
  def _makevars (self):
    import numpy as np
    from fstd2nc.mixins import _axis_type, _var_type
    from collections import OrderedDict
    super(Times,self)._makevars()

    # Annotations for axes
    # Note: 'calendar' attribute is already set to 'noleap' from within
    # fstd2nc.
    time_atts = OrderedDict([('axis','T'),('long_name','time'),('standard_name','time')])

    # Remove degenerate times (single time, encoded as "0").
    for var in self._varlist:
      if 'time' in var.dims:
        ind = var.dims.index('time')
        axis = var.axes[ind]
        if len(axis) == 1 and axis.array[0] == 0:
          var.axes = var.axes[:ind] + var.axes[ind+1:]
          var.record_id = var.record_id.squeeze(axis=ind)

    # Figure out how to interpret the time axis.
    time_type = self._time_type
    # Use heuristic approach to figure out best interpretation if none
    # specicified by the user.
    # 60% of the time, it works every time.
    could_be_from_file = (time_type == 'from-file')
    if time_type == 'auto':
      could_be_year1offset = True
      could_be_yyyymmddhh = True
      could_be_yyyymmdd = True
      could_be_yyyymm = True
      could_be_from_file = True
      for axis in self._iter_axes('time'):
        # Check for consistent increments of values.
        if len(set(np.diff(axis.array))) > 1:
          could_be_year1offset = False
        # Check if value can be interpreted as a date string.
        try:
          np.array(yyyymmddhh_to_date(axis.array))
        except ValueError:
          could_be_yyyymmddhh = False
        try:
          np.array(yyyymmdd_to_date(axis.array))
        except ValueError:
          could_be_yyyymmdd = False
        if len(axis) != 12 or not np.all(axis.array%100 == np.arange(12)+1):
          could_be_yyyymm = False
        # Check if have single time value, and date is in filename.
        if len(axis) > 1 or len(self._files) > 1:
          could_be_from_file = False
        try:
          np.array(filename_to_date(self._files))
        except ValueError:
          could_be_from_file = False

      if could_be_yyyymm: time_type = 'yyyymm'
      elif could_be_yyyymmddhh: time_type = 'yyyymmddhh'
      elif could_be_yyyymmdd: time_type = 'yyyymmdd'
      elif could_be_year1offset: time_type = 'year1-offset'
      else: time_type = 'code'

    # Decode the time values.
    for axis in self._iter_axes('time'):
      if could_be_from_file and axis.array[0] <= 1000:
        axis.array = np.array(filename_to_date(self._files))
      elif time_type == 'year1-offset':
        axis.array = np.array(year1offset_to_date(axis.array))
      elif time_type == 'yyyymmddhh':
        axis.array = np.array(yyyymmddhh_to_date(axis.array))
      elif time_type == 'yyyymmdd':
        axis.array = np.array(yyyymmdd_to_date(axis.array))
      elif time_type == 'yyyymm':
        axis.array = np.array(yyyymm_to_date(axis.array))

      # Annotate the time axis.
      if time_type != 'code':
        axis.atts.update(time_atts)
