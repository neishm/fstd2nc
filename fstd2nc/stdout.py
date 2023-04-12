###############################################################################
# Copyright 2017-2023 - Climate Research Division
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

######################################################################
# Helper functions for printing information / warning / error messages,
# as well as translation support.
######################################################################

# Enable multi-language support.
from gettext import gettext as _
import gettext
from os import path, environ
gettext.bindtextdomain('fstd2nc', path.join(path.dirname(__file__),'locale'))
gettext.textdomain('fstd2nc')
# Check for special CMCLNG environment variable
if environ.get('CMCLNG') == 'francais':
  environ['LANGUAGE'] = 'fr_CA'
del gettext, path, environ

# Module-level variable to control the amount of information printed.
# For Python invocation, only warning and error messages are displayed.
# This is overridden in _fstd2nc_cmdline to include 'info' messages.
streams=('warn','error')
_python = True  # Becomes False for command-line invocation.

from textwrap import TextWrapper
w = TextWrapper(subsequent_indent='    ', break_on_hyphens=False)
del TextWrapper

# Information messages
def info (msg):
  if 'info' in streams:
    print (w.fill(msg))

# How to handle warning messages.
# E.g., can either pass them through warnings.warn, or simply print them.
def warn (msg, _printed=set()):
  from warnings import warn
  if 'warn' not in streams: return
  if _python:
    warn(msg)
  elif msg not in _printed:
    print (_("Warning: %s")%w.fill(msg))
    _printed.add(msg)

# Error messages
def error (msg):
  from sys import exit
  if _python:
    raise Exception(msg)
  elif 'error' in streams:
    print (_("Error: %s")%w.fill(msg))
  exit(1)


