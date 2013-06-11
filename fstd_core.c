#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <string.h>

// Wrap a Fortran call
// Note: change this as needed on your platform.
#define f77name(name) name ## _

extern int c_fnom (int*, char*, char*, int);
extern int c_fclos (int);
extern int c_fstouv (int, char*);
extern int c_fstfrm (int);
extern int c_fstnbr (int);
extern int c_fstinfx (int, int, int*, int*, int*, int, char*, int, int, int, char*, char*);
extern int c_fstprm (int, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, char*, char*, char*, char*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern int c_fstluk (void*, int, int*, int*, int*);
extern int c_fstecr (void*, void*, int, int, int, int, int, int, int, int, int, int, int, char*, char*, char*, char*, int, int, int, int, int, int);
extern int f77name(newdate) (int*, int*, int*, int*);
extern void f77name(convip)(int*, float*, int*, int*, char*, int*);

typedef struct {
  int handle;
  int dateo, deet, npas;
  int ni, nj, nk;
  int nbits, datyp, ip1;
  int ip2, ip3;
  char typvar[2], nomvar[4], etiket[12], grtyp[2];  // grtyp size 2 for better alignment
  int ig1, ig2, ig3;
  int ig4, swa, lng;
  int dltf, ubc, extra1;
  int extra2, extra3;
  PyObject *data_func;
} HEADER;

// Open a file for read access
static PyObject *fstd_open_readonly (PyObject *self, PyObject *args) {
  int iun=0, ier;
  char *filename;
  if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;
  ier = c_fnom (&iun, filename, "STD+RND+R/O", 0);
  if (ier != 0) return NULL;
  c_fstouv (iun, "RND");
  return Py_BuildValue("i", iun);
}

// Open a file for write access
static PyObject *fstd_open_write (PyObject *self, PyObject *args) {
  int iun=0, ier;
  char *filename;
  if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;
  ier = c_fnom (&iun, filename, "STD+RND+R/W", 0);
  if (ier != 0) return NULL;
  c_fstouv (iun, "RND");
  return Py_BuildValue("i", iun);
}

// Close a file
static PyObject *fstd_close (PyObject *self, PyObject *args) {
  int iun, ier;
  if (!PyArg_ParseTuple(args, "i", &iun)) return NULL;
  ier = c_fstfrm(iun);
  if (ier != 0) return NULL;
  ier = c_fclos (iun);
  if (ier != 0) return NULL;
  Py_RETURN_NONE;
}

// Return all record headers in the given file
static PyObject *fstd_get_record_headers (PyObject *self, PyObject *args) {
  int iun, ier, nrec;
  PyObject *header_structure;
  PyArray_Descr *descr;
  PyArrayObject *headers;
  HEADER *h;
  npy_intp dims[1];
  int handle, i;

  if (!PyArg_ParseTuple(args, "i", &iun)) return NULL;
  nrec = c_fstnbr(iun);

  // Allocate the header array
  dims[0] = nrec;
  header_structure = Py_BuildValue("[(s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s)]", "handle", "i4", "dateo", "i4", "deet", "i4", "npas", "i4", "ni", "i4", "nj", "i4", "nk", "i4", "nbits", "i4", "datyp", "i4", "ip1", "i4", "ip2", "i4", "ip3", "i4", "typvar", "a2", "nomvar", "a4", "etiket", "a12", "grtyp", "a2", "ig1", "i4", "ig2", "i4", "ig3", "i4", "ig4", "i4", "swa", "i4", "lng", "i4", "dltf", "i4", "ubc", "i4", "extra1", "i4", "extra2", "i4", "extra3", "i4", "data_func", "O");
  PyArray_DescrConverter (header_structure, &descr);
  Py_DECREF (header_structure);
  headers = (PyArrayObject*) PyArray_SimpleNewFromDescr (1, dims, descr);

  // Populate the header array with data from the file
  h = (HEADER*)headers->data;
  handle = -2;
  for (i = 0; i < nrec; i++) {
    int ni, nj, nk;
    handle = c_fstinfx (handle, iun, &ni, &nj, &nk, -1, "", -1, -1, -1, "", "");
    h->handle = handle;

    memset (h->typvar, ' ', 2);
    memset (h->nomvar, ' ', 4);
    memset (h->etiket, ' ', 12);

    ier = c_fstprm (handle, &h->dateo, &h->deet, &h->npas, &h->ni, &h->nj, &h->nk, &h->nbits, &h->datyp, &h->ip1, &h->ip2, &h->ip3, h->typvar, h->nomvar, h->etiket, h->grtyp, &h->ig1, &h->ig2, &h->ig3, &h->ig4, &h->swa, &h->lng, &h->dltf, &h->ubc, &h->extra1, &h->extra2, &h->extra3);

    if (ier != 0) return NULL;

    h++;
  }

  return (PyObject*)headers;
}

// Read a record from a file
static PyObject *fstd_read_record (PyObject *self, PyObject *args) {
  int ier, handle, ni, nj, nk;
  PyArrayObject *out;
  if (!PyArg_ParseTuple(args, "iO!", &handle, &PyArray_Type, &out)) return NULL;
  ier = c_fstluk (out->data, handle, &ni, &nj, &nk);
  if (ier < 0) return NULL;
  Py_RETURN_NONE;
}

// Write records to a file
static PyObject *fstd_write_records (PyObject *self, PyObject *args) {
  PyObject *field_obj;
  PyArrayObject *headers, *field;
  HEADER *h;
  int iun, nrec, i, ier;
  if (!PyArg_ParseTuple(args, "iO!", &iun, &PyArray_Type, &headers)) return NULL;
  if (PyArray_ITEMSIZE(headers) != sizeof(HEADER)) return NULL;
  printf ("okay\n");
  if (!PyArray_ISCONTIGUOUS(headers)) return NULL;
  nrec = PyArray_SIZE(headers);
  h = (HEADER*)headers->data;

  for (i = 0; i < nrec; i++) {
    if (h->data_func == NULL) return NULL;
    field_obj = PyObject_CallObject (h->data_func, NULL);
    if (field_obj == NULL) return NULL;
    field = (PyArrayObject*)PyArray_EnsureArray(field_obj);
    if (field == NULL) return NULL;
    if (PyArray_SIZE(field) < (h->ni * h->nj * h->nk)) {
      PyErr_SetString (PyExc_TypeError, "Array from data_func is too small");
      Py_DECREF(field);
      return NULL;
    }
    char typvar[3] = {' ',' ','\0'};
    char nomvar[5] = {' ',' ',' ',' ','\0'};
    char etiket[13] = {' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};
    char grtyp[2] = {' ','\0'};
    strncpy (typvar, h->typvar, 2);
    strncpy (nomvar, h->nomvar, 4);
    strncpy (etiket, h->etiket, 12);
    strncpy (grtyp, h->grtyp, 1);
    ier = c_fstecr (field->data, NULL, -h->nbits, iun, h->dateo, h->deet, h->npas, h->ni, h->nj, h->nk, h->ip1, h->ip2, h->ip3, typvar, nomvar, etiket, grtyp, h->ig1, h->ig2, h->ig3, h->ig4, h->datyp, 0);
    if (ier < 0) return NULL;
    Py_DECREF(field);
    h++;
  }
  Py_RETURN_NONE;
}

// Convert CMC timestamps to seconds since 1980-01-01 00:00:00.
static PyObject *stamp2date (PyObject *self, PyObject *args) {
  PyObject *stamp_obj;
  PyArrayObject *date_array, *stamp_array;
  int *date, *stamp, run = 0, mode = 1;
  int ier, i;
  npy_intp n;
  if (!PyArg_ParseTuple(args, "O", &stamp_obj)) return NULL;
  stamp_array = (PyArrayObject*)PyArray_ContiguousFromAny(stamp_obj,NPY_INT,0,0);
  if (stamp_array == NULL) return NULL;
  n = PyArray_SIZE(stamp_array);
  date_array = (PyArrayObject*)PyArray_SimpleNew(1, &n, NPY_INT);
  if (date_array == NULL) {
    Py_DECREF(stamp_array);
    return NULL;
  }
  date = (int*)date_array->data;
  stamp = (int*)stamp_array->data;
  for (i = 0; i < n; i++) {
    ier = f77name(newdate)(date, stamp, &run, &mode);
    if (ier != 0) {
      Py_DECREF(stamp_array);
      Py_DECREF(date_array);
      PyErr_SetString (PyExc_RuntimeError, "Problem encountered with NEWDATE");
      return NULL;
    }
    // Convert from 5-second intervals
    *date *= 5;
    date++; stamp++;
  }
  Py_DECREF(stamp_array);
  return (PyObject*)date_array;
}

// Convert seconds since 1980-01-01 00:00:00 to CMC timestamps.
static PyObject *date2stamp (PyObject *self, PyObject *args) {
  PyObject *date_obj;
  PyArrayObject *date_array, *stamp_array;
  int *date, *stamp, run = 0, mode = -1;
  int ier, i;
  npy_intp n;
  if (!PyArg_ParseTuple(args, "O", &date_obj)) return NULL;
  date_array = (PyArrayObject*)PyArray_ContiguousFromAny(date_obj,NPY_INT,0,0);
  if (date_array == NULL) return NULL;
  n = PyArray_SIZE(date_array);
  stamp_array = (PyArrayObject*)PyArray_SimpleNew(1, &n, NPY_INT);
  if (stamp_array == NULL) {
    Py_DECREF(date_array);
    return NULL;
  }
  date = (int*)date_array->data;
  stamp = (int*)stamp_array->data;
  for (i = 0; i < n; i++) {
    // Convert from 5-second intervals
    int current_date = *date / 5;
    ier = f77name(newdate)(&current_date, stamp, &run, &mode);
    if (ier != 0) {
      Py_DECREF(stamp_array);
      Py_DECREF(date_array);
      PyErr_SetString (PyExc_RuntimeError, "Problem encountered with NEWDATE");
      return NULL;
    }
    date++; stamp++;
  }
  Py_DECREF(date_array);
  return (PyObject*)stamp_array;
}

// Decode vertical levels
static PyObject *decode_levels (PyObject *self, PyObject *args) {
  PyObject *ip1_obj;
  PyArrayObject *ip1_array, *z_array;
  float *z;
  int *ip1, kind, mode = -1, flag = 0;
  int i;
  npy_intp n;
  if (!PyArg_ParseTuple(args, "O", &ip1_obj)) return NULL;
  ip1_array = (PyArrayObject*)PyArray_ContiguousFromAny(ip1_obj,NPY_INT,0,0);
  if (ip1_array == NULL) return NULL;
  n = PyArray_SIZE(ip1_array);
  z_array = (PyArrayObject*)PyArray_SimpleNew(1, &n, NPY_FLOAT32);
  if (z_array == NULL) {
    Py_DECREF(ip1_array);
    return NULL;
  }
  ip1 = (int*)ip1_array->data;
  z = (float*)z_array->data;
  for (i = 0; i < n; i++) {
    f77name(convip)(ip1++, z++, &kind, &mode, "", &flag);
  }
  PyObject *ret = Py_BuildValue("(O,i)", z_array, kind);
  Py_DECREF (ip1_array);
  Py_DECREF (z_array);
  return ret;
}

// Encode vertical levels
static PyObject *encode_levels (PyObject *self, PyObject *args) {
  PyObject *z_obj;
  PyArrayObject *ip1_array, *z_array;
  float *z;
  int *ip1, kind, mode = 2, flag = 0;
  int i;
  npy_intp n;
  if (!PyArg_ParseTuple(args, "Oi", &z_obj, &kind)) return NULL;
  z_array = (PyArrayObject*)PyArray_ContiguousFromAny(z_obj,NPY_FLOAT32,0,0);
  if (z_array == NULL) return NULL;
  n = PyArray_SIZE(z_array);
  ip1_array = (PyArrayObject*)PyArray_SimpleNew(1, &n, NPY_INT);
  if (ip1_array == NULL) {
    Py_DECREF(z_array);
    return NULL;
  }
  ip1 = (int*)ip1_array->data;
  z = (float*)z_array->data;
  for (i = 0; i < n; i++) {
    f77name(convip)(ip1++, z++, &kind, &mode, "", &flag);
  }
  Py_DECREF (z_array);
  return (PyObject*)ip1_array;
}

static PyMethodDef FST_Methods[] = {
  {"open_readonly", fstd_open_readonly, METH_VARARGS, "Open an FSTD file for read access"},
  {"open_write", fstd_open_write, METH_VARARGS, "Open an FSTD file for write access"},
  {"close", fstd_close, METH_VARARGS, "Close an FSTD file"},
  {"get_record_headers", fstd_get_record_headers, METH_VARARGS, "Get all record headers from an FSTD file"},
  {"read_record", fstd_read_record, METH_VARARGS, "Read a record into the given numpy array"},
  {"write_records", fstd_write_records, METH_VARARGS, "Write a set of records into a given FSTD file"},
  {"stamp2date", stamp2date, METH_VARARGS, "Convert CMC timestamps to seconds since 1980-01-01 00:00:00"},
  {"date2stamp", date2stamp, METH_VARARGS, "Convert seconds since 1980-01-01 00:00:00 to a CMC timestamp"},
  {"decode_levels", decode_levels, METH_VARARGS, "Decode vertical levels"},
  {"encode_levels", encode_levels, METH_VARARGS, "Encode vertical levels"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initfstd_core(void) {
  (void) Py_InitModule("fstd_core", FST_Methods);
  import_array();
}

