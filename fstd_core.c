#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <string.h>

extern int c_fnom (int*, char*, char*, int);
extern int c_fclos (int);
extern int c_fstouv (int, char*);
extern int c_fstfrm (int);
extern int c_fstnbr (int);
extern int c_fstinfx (int, int, int*, int*, int*, int, char*, int, int, int, char*, char*);
extern int c_fstprm (int, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, char*, char*, char*, char*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern int c_fstluk (void*, int, int*, int*, int*);

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

static PyMethodDef FST_Methods[] = {
  {"open_readonly", fstd_open_readonly, METH_VARARGS, "Open an FSTD file for read access"},
  {"close", fstd_close, METH_VARARGS, "Close an FSTD file"},
  {"get_record_headers", fstd_get_record_headers, METH_VARARGS, "Get all record headers from an FSTD file"},
  {"read_record", fstd_read_record, METH_VARARGS, "Read a record into the given numpy array"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initfstd_core(void) {
  (void) Py_InitModule("fstd_core", FST_Methods);
  import_array();
}

