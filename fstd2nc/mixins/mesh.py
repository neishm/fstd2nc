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

from fstd2nc.stdout import _, info, warn, error
from fstd2nc.mixins import BufferBase

#################################################
# Mixin for handling triangular mesh grids.
#
# https://ugrid-conventions.github.io/ugrid-conventions/
#

class Mesh (BufferBase):
  def __init__ (self, *args, **kwargs):
    # Tell the decoder not to process horizontal records as variables.
    #self._maybe_meta_records = self._maybe_meta_records + (b'##',)
    super(Mesh,self).__init__(*args,**kwargs)

  # Add horizontal coordinate info to the data.
  def _makevars (self):
    from fstd2nc.mixins import _iter_type, _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    import numpy as np

    super(Mesh,self)._makevars()

    # Scan for mesh coordinates.
    three = _dim_type ('Three', 3)
    mesh_indices = OrderedDict()
    lats = OrderedDict()
    lons = OrderedDict()
    for var in self._varlist:
      # Look for face node indices ('##' field).
      if var.name == '##':
        var.name = 'mesh_face_nodes'
        # Reshape to (nFaces,3)
        if 'i' in var.dims:
          ind = var.dims.index('i')
          faces = _dim_type ('nFaces', var.shape[ind]//3)
          if var.shape[ind] % 3 == 0:
            var.axes = var.axes[:ind] + [faces, three] + var.axes[ind+1:]
            shape = var.record_id.shape
        # Remove degenerate height dimension.
        if 'height' in var.dims:
          ind = var.dims.index('height')
          var.axes = var.axes[:ind] + var.axes[ind+1:]
          var.record_id = var.record_id.reshape(shape[:ind]+shape[ind+1:])
        mesh_indices[(var.atts['ip1'],var.atts['ip2'])] = var
        # Remove coordinate arrays.
        var.atts.pop('coordinates',None)
      # Look for lat/lon coordinates.
      if var.atts.get('grtyp',None) == 'M':
        for coord in var.atts.get('coordinates',[]):
          if coord.name == 'lat':
            lats[var.atts['ig1'],var.atts['ig2']] = coord
          if coord.name == 'lon':
            lons[var.atts['ig1'],var.atts['ig2']] = coord

    if len(mesh_indices) == 0: return

    # Create toplogy variable(s)
    topos = OrderedDict()
    for key in mesh_indices.keys():
      coordinates = []
      if key in lats: coordinates.append(lats[key])
      if key in lons: coordinates.append(lons[key])
      atts = OrderedDict()
      atts['cf_role'] = 'mesh_topology'
      atts['topology_dimension'] = 2
      atts['node_coordinates'] = coordinates
      atts['face_node_connectivity'] = mesh_indices[key]
      topo = _var_type ('Mesh', atts, (), np.int32(0))
      topos[key] = topo

    # Add toplogy information to data variables.
    # Also, update 'i' dimension to be nNodes for the data variables.
    reftimes = OrderedDict()
    times = OrderedDict()
    for var in self._varlist:
      if var.atts.get('grtyp',None) == 'M':
        key = (var.atts['ig1'], var.atts['ig2'])
        if key in topos:
          var.atts['mesh'] = topos[key]
          var.atts['location'] = 'node'
        if 'i' in var.dims:
          ind = var.dims.index('i')
          var.axes[ind].name = 'nNodes'
        # Get a copy of the forecast times (to adjust the coordinates).
        if key not in reftimes:
          for coord in var.atts.get('coordinates',[]):
            if coord.name == 'reftime':
              reftimes[key] = coord
          if 'time' in var.dims:
            times[key] = var.axes[var.dims.index('time')]

    # Need to fix the time axis for the lat / lon / triangle indices
    # (from '^^', '>>', '##' records).
    # They may not include the forecast period.
    for key in topos.keys():
      if key in times and key in reftimes:
        for category in mesh_indices, lats, lons:
          if key in category:
            var = category[key]
            if 'time' in var.dims:
              ind = var.dims.index('time')
              if np.all(var.axes[ind].array == reftimes[key].array):
                var.axes[ind] = times[key]

