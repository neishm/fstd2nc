import fstd2nc
from fstd2nc import extra as fstd2nc_extra # to avoid direct extra.py
from fstd2nc.mixins import _iter_type, _axis_type, _dim_type # For type checking
import numpy

class FSTDInspector:
    def __init__(self, filepaths, **buffer_kwargs):
        self.buffer = fstd2nc.Buffer(filepaths, **buffer_kwargs)
        self.buffer._makevars()

    def get_variables_info(self):
        variables_metadata = []
        for var in self.buffer._varlist: # or self.buffer._iter_objects()
            if isinstance(var, _iter_type):
                name = var.name
                attributes = var.atts.copy()
                dtype_str = var.dtype.str
                dimensions = var.dims
                shape = var.shape
                
                coordinates_info = []
                for axis_obj in var.axes:
                    if isinstance(axis_obj, _axis_type):
                        coordinates_info.append({
                            'name': axis_obj.name,
                            'data': axis_obj.array,
                            'attributes': axis_obj.atts.copy()
                        })
                    elif isinstance(axis_obj, _dim_type):
                        coordinates_info.append({
                            'name': axis_obj.name,
                            'length': axis_obj.size # Assuming _dim_type has a 'size' or 'length' attribute
                        })
                
                variables_metadata.append({
                    'name': name,
                    'attributes': attributes,
                    'dtype': dtype_str,
                    'dimensions': dimensions,
                    'shape': shape,
                    'coordinates': coordinates_info
                })
        return variables_metadata

    def get_chunk_info(self, variable_name):
        var = None
        for v in self.buffer._varlist:
            if v.name == variable_name and isinstance(v, _iter_type):
                var = v
                break
        
        if var is None:
            raise ValueError(f"Variable '{variable_name}' not found.")

        chunks_data = []
        for idx, rec_idx in numpy.ndenumerate(var.record_id):
            if rec_idx < 0:  # No data for this logical chunk
                continue
            
            file_path = self.buffer._files[self.buffer._headers['file_id'][rec_idx]]
            offset = self.buffer._headers['address'][rec_idx]
            length = self.buffer._headers['length'][rec_idx]
            logical_coords = idx  # N-D index from ndenumerate
            
            chunks_data.append({
                'file_path': file_path,
                'offset': offset,
                'length': length,
                'logical_coords': logical_coords
            })
        return chunks_data

    def decode_chunk(self, file_path, offset, length):
        with open(file_path, 'rb') as f:
            f.seek(offset)
            raw_bytes_segment = f.read(length)
        
        # Ensure raw_bytes_segment is a NumPy array of dtype 'B'
        raw_bytes_segment_np = numpy.frombuffer(raw_bytes_segment, dtype='B')
        
        decoded_data = fstd2nc_extra.decode(raw_bytes_segment_np)
        return decoded_data
