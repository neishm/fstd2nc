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

from fstd2nc.stdout import _, info, warn, error
from fstd2nc.mixins import BufferBase

#################################################
# Mixin for handling triangular mesh grids.
#
# https://ugrid-conventions.github.io/ugrid-conventions/
#

class Mesh (BufferBase):
  # Tell the decoder not to process horizontal records as variables.
  @classmethod
  def _maybe_meta_records(cls):
    return super(Mesh,cls)._maybe_meta_records() + (b'##',)

  def __init__ (self, *args, **kwargs):
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
      if var.atts.get('grtyp',None) != 'M': continue
      key = (var.atts['ig1'], var.atts['ig2'])
      if key not in mesh_indices:
        # Look for face node indices ('##' field).
        mesh = self._find_coord(var,b'##  ')  # Borrowed from xycoords.
        d = mesh['d'].transpose().squeeze()
        if (d.ndim == 1):
          faces = _dim_type ('nFaces', len(d)//3)
          d = d.reshape(len(d)//3,-1)
          atts = OrderedDict([("long_name","Vertices index of triangular meshe (TIN)"),("description","Index of the vertices position encoded in the ^^ >> fields which constitutes a triangular mesh")])
          mesh = _var_type('mesh_face_nodes',atts, [faces,three], d)
          mesh_indices[key] = mesh
      # Look for lat/lon coordinates.
      for coord in var.atts.get('coordinates',[]):
        if coord.name == 'lat':
          lats[key] = coord
        if coord.name == 'lon':
          lons[key] = coord

    if len(mesh_indices) == 0: return

    # Create toplogy variable(s)
    topos = OrderedDict()
    for key in mesh_indices.keys():
      coordinates = []
      if key in lats: coordinates.append(lats[key])
      if key in lons: coordinates.append(lons[key])
      atts = OrderedDict()
      atts['cf_role'] = 'mesh_topology'
      atts['topology_dimension'] = np.int32(2)
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

