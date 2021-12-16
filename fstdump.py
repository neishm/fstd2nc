###############################################################################
# Copyright 2017-2021 - Climate Research Division
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


from __future__ import print_function
from fstd2nc.stdout import _, info, warn, error
import fstd2nc
from fstd2nc import Buffer

#################################################
# Dump the metadata for an RPN standard file to stdout (without conversion
# step).  Displays in format similar to "ncdump" for netCDF files.

def fstdump (buffer_type=Buffer):
  from argparse import ArgumentParser
  from sys import stdout, argv
  from os.path import exists
  from rpnpy.librmn.fstd98 import FSTDError, fstopt
  parser = ArgumentParser(description=_("Display RPN standard file (FSTD) metadata in a structured format."))
  parser.add_argument('infile', metavar='<fstd_file>', help=_('The RPN standard file to query.'))
  parser.add_argument('-v', nargs='*', metavar='NOMVAR,...', help=_('Display the values for the specified variable.'))
  buffer_type._cmdline_args(parser)
  args = parser.parse_args()
  buffer_type._check_args(parser, args)
  args = vars(args)
  infile = args.pop('infile')
  nomvars = args.pop('v')

  # Apply message level criteria.
  fstopt ('MSGLVL',6)
  fstd2nc.stdout.streams = ('error',)

  try:
    buf = buffer_type(infile, **args)
  except FSTDError:
    error (_("problem opening one or more input files."))

  # Get the metadata in a netCDF-like structure.
  x = buf.to_xarray()

  # Print header info.
  x.info()
  print ()

  # Print variable data.
  if nomvars is not None:
    nomvars = [n for nomvar in nomvars for n in nomvar.split(',')]
    for nomvar in nomvars:
      if nomvar in x:
        values = x[nomvar].values
        if 'bytes' in values.dtype.name:
          values = values.astype(str)
        print ("\n%s = %s ;"%(nomvar,values))

if __name__ == '__main__':
  fstdump()

