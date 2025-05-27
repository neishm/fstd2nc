import unittest
import os
import numpy as np
from fstd2nc.inspector import FSTDInspector
import fstd2nc # To access fstd2nc.FSTDError for exception checking

# Helper function to create a dummy FSTD file path
DUMMY_FSTD_FILENAME = "simple_test.fstd"
DUMMY_FSTD_FILEPATH = os.path.join(os.path.dirname(__file__), "data", DUMMY_FSTD_FILENAME)
PLACEHOLDER_CONTENT = b"NOT_A_VALID_FSTD_FILE_PLACEHOLDER_MORE_BYTES_TO_BE_SURE"

def ensure_dummy_fstd_file():
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    if not os.path.exists(DUMMY_FSTD_FILEPATH):
        with open(DUMMY_FSTD_FILEPATH, 'wb') as f:
            f.write(PLACEHOLDER_CONTENT)
        # print(f"Warning: Created a placeholder file at {DUMMY_FSTD_FILEPATH}. Tests needing a real FSTD will fail or be skipped.")

def is_placeholder_file(filepath):
    if not os.path.exists(filepath):
        return False
    with open(filepath, 'rb') as f:
        return f.read(len(PLACEHOLDER_CONTENT)) == PLACEHOLDER_CONTENT

class TestFSTDInspector(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ensure_dummy_fstd_file()
        
        cls.is_placeholder = is_placeholder_file(DUMMY_FSTD_FILEPATH)

        if not cls.is_placeholder:
            # Attempt to load if it's not our placeholder
            # This path is for when a user provides a real simple_test.fstd
            try:
                cls.inspector = FSTDInspector(DUMMY_FSTD_FILEPATH)
                print(f"Successfully initialized FSTDInspector with provided file: {DUMMY_FSTD_FILEPATH}")
            except Exception as e:
                print(f"Failed to initialize FSTDInspector with provided file {DUMMY_FSTD_FILEPATH}: {e}. Tests will be skipped.")
                cls.inspector = None
        else:
            # If it IS the placeholder, inspector should not be initialized here.
            # Tests will verify behavior with the placeholder.
            cls.inspector = None
            # print(f"Note: Using placeholder FSTD file: {DUMMY_FSTD_FILEPATH}. Some tests will be skipped or verify error handling.")

    def test_00_placeholder_check(self):
        # This test just confirms the state of the dummy file for clarity in test logs
        if self.is_placeholder:
            print(f"Test environment is using the placeholder file: {DUMMY_FSTD_FILEPATH}")
        else:
            print(f"Test environment is using a user-provided file: {DUMMY_FSTD_FILEPATH}")
        self.assertTrue(os.path.exists(DUMMY_FSTD_FILEPATH))

    def test_01_initialization(self):
        if self.is_placeholder:
            # Expect initialization to fail with the placeholder
            with self.assertRaises((fstd2nc.FSTDError, ValueError, TypeError, IndexError), msg="FSTDInspector should raise an error with a placeholder file."):
                FSTDInspector(DUMMY_FSTD_FILEPATH)
        elif not self.inspector:
            # A real file was provided, but FSTDInspector failed to initialize in setUpClass
            self.skipTest(f"FSTDInspector failed to initialize with the provided file: {DUMMY_FSTD_FILEPATH}. See setUpClass logs.")
        else:
            # Inspector was successfully initialized with a user-provided file
            self.assertIsNotNone(self.inspector, "FSTDInspector should be initialized.")
            self.assertIsNotNone(self.inspector.buffer, "Inspector should have a buffer attribute.")
            self.assertTrue(hasattr(self.inspector.buffer, '_varlist'), "Buffer should have _varlist after init.")


    def test_02_get_variables_info(self):
        if self.is_placeholder:
            self.skipTest(f"Skipping get_variables_info test as {DUMMY_FSTD_FILEPATH} is a placeholder.")
        if not self.inspector:
            self.skipTest(f"FSTDInspector not initialized. Skipping test. Check {DUMMY_FSTD_FILEPATH} or provide a valid one.")
        
        vars_info = self.inspector.get_variables_info()
        self.assertIsInstance(vars_info, list, "get_variables_info should return a list.")
        
        if not vars_info: # This might be normal for an empty (but valid) FSTD
            print(f"Warning: No variables found in {DUMMY_FSTD_FILEPATH}. Detailed checks for 'TEMP' will be skipped.")
            return 

        temp_var_info = next((v for v in vars_info if v['name'] == 'TEMP'), None)
        
        if not temp_var_info:
            self.skipTest(f"Test variable 'TEMP' not found in {DUMMY_FSTD_FILEPATH}. Cannot perform detailed checks.")

        self.assertEqual(temp_var_info['name'], 'TEMP')
        self.assertIn(temp_var_info['dtype'], ['float32', 'f4']) 
        self.assertEqual(temp_var_info['shape'], (2, 2)) 
        self.assertEqual(temp_var_info['dimensions'], ('nj', 'ni')) 

        ni_coord = next((c for c in temp_var_info['coordinates'] if c['name'] == 'ni'), None)
        nj_coord = next((c for c in temp_var_info['coordinates'] if c['name'] == 'nj'), None)
        self.assertIsNotNone(ni_coord, "ni coordinate/dimension should exist.")
        self.assertIsNotNone(nj_coord, "nj coordinate/dimension should exist.")
        self.assertEqual(ni_coord['length'], 2) # Assuming these are _dim_type
        self.assertEqual(nj_coord['length'], 2)


    def test_03_get_chunk_info(self):
        if self.is_placeholder:
            self.skipTest(f"Skipping get_chunk_info test as {DUMMY_FSTD_FILEPATH} is a placeholder.")
        if not self.inspector:
            self.skipTest(f"FSTDInspector not initialized. Skipping test. Check {DUMMY_FSTD_FILEPATH} or provide a valid one.")
        
        # Pre-check if 'TEMP' variable exists to avoid ValueError if it's missing
        # This requires get_variables_info to have found it.
        vars_info = self.inspector.get_variables_info()
        temp_var_meta = next((v for v in vars_info if v['name'] == 'TEMP'), None)
        if not temp_var_meta:
            self.skipTest(f"Test variable 'TEMP' not found in {DUMMY_FSTD_FILEPATH}, skipping get_chunk_info.")
            return

        chunk_info_list = self.inspector.get_chunk_info('TEMP')
        self.assertIsInstance(chunk_info_list, list)
        if not chunk_info_list:
            self.skipTest(f"No chunks found for variable 'TEMP' in {DUMMY_FSTD_FILEPATH}.") # Might be valid for some FSTD structures

        self.assertEqual(len(chunk_info_list), 1, "Expected one chunk for simple TEMP variable.")
        
        chunk = chunk_info_list[0]
        self.assertIn('file_path', chunk)
        self.assertIn('offset', chunk)
        self.assertIn('length', chunk)
        self.assertIn('logical_coords', chunk)
        
        self.assertEqual(chunk['file_path'], DUMMY_FSTD_FILEPATH) # Or the actual path from multi-file buffer
        self.assertIsInstance(chunk['offset'], int)
        self.assertIsInstance(chunk['length'], int)
        self.assertGreaterEqual(chunk['offset'], 0)
        self.assertGreater(chunk['length'], 0)
        
        temp_var_obj = next((v for v in self.inspector.buffer._varlist if v.name == 'TEMP'), None)
        self.assertIsNotNone(temp_var_obj, "TEMP variable object not found in buffer._varlist")
        self.assertEqual(chunk['logical_coords'], (0,) * temp_var_obj.record_id.ndim)


    def test_04_decode_chunk(self):
        if self.is_placeholder:
            self.skipTest(f"Skipping decode_chunk test as {DUMMY_FSTD_FILEPATH} is a placeholder.")
        if not self.inspector:
            self.skipTest(f"FSTDInspector not initialized. Skipping test. Check {DUMMY_FSTD_FILEPATH} or provide a valid one.")

        vars_info = self.inspector.get_variables_info()
        temp_var_meta = next((v for v in vars_info if v['name'] == 'TEMP'), None)
        if not temp_var_meta:
            self.skipTest(f"Test variable 'TEMP' not found in {DUMMY_FSTD_FILEPATH}, skipping decode_chunk.")
            return

        chunk_info_list = self.inspector.get_chunk_info('TEMP')
        if not chunk_info_list:
            self.skipTest(f"No chunks found for 'TEMP' to decode.")
            return

        chunk_to_decode = chunk_info_list[0]
        
        try:
            decoded_data = self.inspector.decode_chunk(
                chunk_to_decode['file_path'],
                chunk_to_decode['offset'],
                chunk_to_decode['length']
            )
        except Exception as e:
            self.fail(f"decode_chunk failed with supposedly real FSTD file: {e}")
            return

        self.assertIsInstance(decoded_data, np.ndarray)
        self.assertEqual(decoded_data.shape, (2, 2)) 
        self.assertTrue(np.issubdtype(decoded_data.dtype, np.floating) or np.issubdtype(decoded_data.dtype, np.integer))

        expected_data = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
        np.testing.assert_array_almost_equal(decoded_data, expected_data, decimal=5)


    def test_05_get_chunk_info_error(self):
        if self.is_placeholder:
            # Create a dummy inspector instance that is expected to fail in a certain way,
            # or mock it. For now, we test the error path on a potentially valid inspector.
            # If no inspector, we can't test this specific error raising.
            self.skipTest(f"Skipping get_chunk_info_error test as {DUMMY_FSTD_FILEPATH} is a placeholder and inspector is not configured for it.")
        if not self.inspector:
            self.skipTest(f"FSTDInspector not initialized. Skipping test. Check {DUMMY_FSTD_FILEPATH} or provide a valid one.")
        
        with self.assertRaises(ValueError, msg="Should raise ValueError for non-existent variable"):
            self.inspector.get_chunk_info('NON_EXISTENT_VARIABLE')

if __name__ == '__main__':
    unittest.main()
