#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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
extern int f77name(ig_to_hybref)(int*, int*, int*, int*, float*, float*, float*, float*);
extern int f77name(hybref_to_ig)(int*, int*, int*, int*, float*, float*, float*, float*);
extern void f77name(cxgaig)(char*, int*, int*, int*, int*, float*, float*, float*, float*);
extern void f77name(cigaxg)(char*, float*, float*, float*, float*, int*, int*, int*, int*);
extern int c_ezqkdef(int ni, int nj, char *grtyp, 
                 int ig1, int ig2, int ig3, int ig4, int iunit);
extern int c_ezgdef_fmem(int ni, int nj, char *grtyp, char *grref, 
                int ig1ref, int ig2ref, int ig3ref, int ig4ref, float *ax, float *ay);
extern int c_gdll (int gdid,  float *lat, float *lon);
extern int c_gdgaxes(int gdid,  float *ax, float *ay);
extern int c_gdrls(int gdid);

// C struct for accessing records
typedef struct {
  int pad;
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

// Corresponding numpy descr (to be initialized)
static PyArray_Descr *descr;


// Data type for holding an FSTD unit.
// Allows the file to be closed when all references are gone.
typedef struct {
  PyObject_HEAD
  int iun;
} FSTD_Unit_Object;

static void FSTD_Unit_dealloc (FSTD_Unit_Object *self) {
  c_fstfrm(self->iun);
  c_fclos (self->iun);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject FSTD_Unit_Type = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "fstd_core.FSTD_Unit",     /*tp_name*/
  sizeof(FSTD_Unit_Object),  /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)FSTD_Unit_dealloc, /*tp_dealloc*/
//  0,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,        /*tp_flags*/
  "FSTD file unit",          /* tp_doc */
};



// Python interface for reading a single record
// (function closure for encapsulating the file and record handle).
typedef struct {
  PyObject_HEAD
  FSTD_Unit_Object *file;
  int handle;
  npy_intp dims[3];  // Dimensions of the record
  int typenum;
} RecordGetter_Object;

static PyObject *RecordGetter_call (PyObject *self, PyObject *args, PyObject *kwargs) {
  RecordGetter_Object *o = (RecordGetter_Object*)self;
  PyArrayObject *out = (PyArrayObject*)PyArray_SimpleNew(3, o->dims, o->typenum);
  if (out == NULL) return NULL;
  int ni, nj, nk;
  c_fstluk (out->data, o->handle, &ni, &nj, &nk);
  return (PyObject*)out;
}

static void RecordGetter_dealloc (RecordGetter_Object *self) {
  Py_DECREF(self->file);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject RecordGetter_Type = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "fstd_core.RecordGetter",  /*tp_name*/
  sizeof(RecordGetter_Object), /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)RecordGetter_dealloc,      /*tp_dealloc*/
//  0,                         /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  RecordGetter_call,         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,        /*tp_flags*/
  "Closures for reading particular FSTD records",           /* tp_doc */
};


// Return all record headers in the given file.
// Include hooks to data getters.
static PyObject *fstd_read_records (PyObject *self, PyObject *args) {
  char *filename;
  int iun = 0, ier, nrec;

  PyArrayObject *headers;
  HEADER *h;
  npy_intp dims[1];
  int handle, i;

  if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;
  ier = c_fnom (&iun, filename, "STD+RND+R/O", 0);
  if (ier != 0) return NULL;
  nrec = c_fstouv (iun, "RND");
  if (nrec < 0) return NULL;
  FSTD_Unit_Object *file = (FSTD_Unit_Object*) PyType_GenericNew (&FSTD_Unit_Type, NULL, NULL);
  if (file == NULL) return NULL;
  file->iun = iun;

  // Allocate the header array
  dims[0] = nrec;
  Py_INCREF(descr);
  headers = (PyArrayObject*) PyArray_SimpleNewFromDescr (1, dims, descr);
  if (headers == NULL) return NULL;

  // Populate the header array with data from the file
//  printf ("sizeof(HEADER) = %d; PyArray_ITEMSIZE(headers) = %d\n", sizeof(HEADER), PyArray_ITEMSIZE(headers));
  if (sizeof(HEADER) != PyArray_ITEMSIZE(headers)) return NULL;
  h = (HEADER*)headers->data;
  handle = -2;
  for (i = 0; i < nrec; i++) {
    int ni, nj, nk;
    handle = c_fstinfx (handle, iun, &ni, &nj, &nk, -1, "", -1, -1, -1, "", "");
    if (handle < 0) return NULL;

    memset (h->typvar, ' ', 2);
    memset (h->nomvar, ' ', 4);
    memset (h->etiket, ' ', 12);

    ier = c_fstprm (handle, &h->dateo, &h->deet, &h->npas, &h->ni, &h->nj, &h->nk, &h->nbits, &h->datyp, &h->ip1, &h->ip2, &h->ip3, h->typvar, h->nomvar, h->etiket, h->grtyp, &h->ig1, &h->ig2, &h->ig3, &h->ig4, &h->swa, &h->lng, &h->dltf, &h->ubc, &h->extra1, &h->extra2, &h->extra3);

    if (ier != 0) return NULL;

    RecordGetter_Object *func = (RecordGetter_Object*) PyType_GenericNew (&RecordGetter_Type, NULL, NULL);
    if (func == NULL) return NULL;
    Py_INCREF(file);
    func->file = file;
    func->handle = handle;
    func->dims[0] = nk;
    func->dims[1] = nj;
    func->dims[2] = ni;
//dtypes = {1:'float', 2:'uint', 3:'a', 4:'int', 5:'float', 134:'float', 130:'uint', 132:'int', 133:'float'}
    func->typenum = -1;
    switch (h->datyp) {
      case 1:
      case 5:
      case 134:
      case 133:
        func->typenum = h->nbits > 32 ? NPY_FLOAT64 : NPY_FLOAT32;
        break;
      case 2:
      case 130:
        func->typenum = h->nbits > 32 ? NPY_UINT64 : NPY_UINT32;
        break;
      case 3:
        func->typenum = NPY_UINT8;
        break;
      case 4:
      case 132:
        func->typenum = h->nbits > 32 ? NPY_INT64 : NPY_INT32;
        break;
      default:
        return NULL;
    }
    h->data_func = (PyObject*)func;

    h++;
  }

  Py_DECREF (file);  // File is referenced in all record getters.

//  Py_INCREF(headers);  // Force no cleanup
  return (PyObject*)headers;
}


// Write records to a file
static PyObject *fstd_write_records (PyObject *self, PyObject *args) {
  char *filename;
  PyObject *field_obj;
  PyArrayObject *headers, *field;
  HEADER *h;
  int iun = 0, nrec, i, ier;
  if (!PyArg_ParseTuple(args, "sO!", &filename, &PyArray_Type, &headers)) return NULL;
  ier = c_fnom (&iun, filename, "STD+RND+R/W", 0);
  if (ier != 0) return NULL;
  c_fstouv (iun, "RND");

  // Ensure the record array is contiguous (and of the right type)
  if (headers->descr != descr) {
    PyErr_SetString (PyExc_ValueError, "Invalid record format");
    return NULL;
  }
  headers = PyArray_GETCONTIGUOUS(headers);
  if (headers == NULL) return NULL;

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

  Py_DECREF(headers);

  c_fstfrm (iun);
  c_fclos (iun);

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
    ier = 0;
    if (*stamp == 0) *date = 0;
    else ier = f77name(newdate)(date, stamp, &run, &mode);
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
  long long *date;
  int *stamp, run = 0, mode = -1;
  int ier, i;
  npy_intp n;
  if (!PyArg_ParseTuple(args, "O", &date_obj)) return NULL;
  date_array = (PyArrayObject*)PyArray_ContiguousFromAny(date_obj,NPY_LONGLONG,0,0);
  if (date_array == NULL) return NULL;
  n = PyArray_SIZE(date_array);
  stamp_array = (PyArrayObject*)PyArray_SimpleNew(1, &n, NPY_INT);
  if (stamp_array == NULL) {
    Py_DECREF(date_array);
    return NULL;
  }
  date = (long long*)date_array->data;
  stamp = (int*)stamp_array->data;
  for (i = 0; i < n; i++) {
    ier = 0;
    // Convert from 5-second intervals
    int current_date = *date / 5;
    if (current_date == 0) *stamp = 0;
    else ier = f77name(newdate)(&current_date, stamp, &run, &mode);
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

// Get hybrid A and B coefficients
static PyObject *get_hybrid_a_b (PyObject *self, PyObject *args) {
  PyArrayObject *hy_record;
  int ig1, ig2, ig3, ig4;
  float ptop, rcoef, pref, x1, x2;
  PyObject *levels_obj;
  PyArrayObject *levels;
  PyArrayObject *A, *B;
  npy_intp n;
  int i;

  // Parse inputs
  if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &hy_record, &levels_obj)) return NULL;
  if (hy_record->descr != descr) return NULL;
  hy_record = (PyArrayObject*) PyArray_GETCONTIGUOUS(hy_record);
  levels = (PyArrayObject*)PyArray_ContiguousFromAny(levels_obj,NPY_FLOAT32,1,1);
  if (levels == NULL) return NULL;

  // Get HY parameters
  HEADER *hy = (HEADER*)hy_record->data;
  ig1 = hy->ig1;
  ig2 = hy->ig2;
  ig3 = hy->ig3;
  ig4 = hy->ig4;

  // Get ptop
  {
    int kind, mode = -1, flag = 0;
    f77name(convip)(&(hy->ip1), &ptop, &kind, &mode, "", &flag);
  }

  // Get rcoef and pref
  f77name(ig_to_hybref)(&ig1, &ig2, &ig3, &ig4, &rcoef, &pref, &x1, &x2);

  // Generate A and B
  n = PyArray_SIZE(levels);
  A = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_FLOAT32);
  B = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_FLOAT32);
  float etatop = ptop / pref;
  for (i = 0; i < n; i++) {
    float *a = (float*)(A->data) + i;
    float *b = (float*)(B->data) + i;
    float eta = ((float*)levels->data)[i];
    *b = pow((eta - etatop) / (1 - etatop), rcoef);
    *a = pref * (eta - *b);
  }

  PyObject *ref = Py_BuildValue("(f,f,f,O,O)", ptop, rcoef, pref, A, B);
  Py_DECREF(A);
  Py_DECREF(B);
  Py_DECREF(hy_record);
  Py_DECREF(levels);
  return ref;
}

// Get zeta (log hybrid) table info (full A and B on momentum & thermodynamic levels)
static PyObject *get_loghybrid_table (PyObject *self, PyObject *args) {
  PyArrayObject *bangbang_record;
  PyObject *table_obj;
  PyArrayObject *table_array;
  double *table;
  PyArrayObject *IP1_m, *A_m, *B_m, *IP1_t, *A_t, *B_t;
  int *ip1_m, *ip1_t;
  double *a_m, *b_m, *a_t, *b_t;
  npy_intp n_m, n_t;
  int i;
  int kind, version;
  double ptop, pref, rcoef1, rcoef2;
  unsigned char ref_name[4];

  // Parse inputs
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &bangbang_record)) return NULL;
  if (bangbang_record->descr != descr) return NULL;
  bangbang_record = (PyArrayObject*) PyArray_GETCONTIGUOUS(bangbang_record);

  // Read !! table info
  table_obj = PyObject_CallObject(((HEADER*)PyArray_DATA(bangbang_record))[0].data_func,NULL);
  Py_DECREF(bangbang_record);
  if (table_obj == NULL) return NULL;
  table_array = (PyArrayObject*)PyArray_ContiguousFromAny(table_obj,NPY_FLOAT64,3,3);
  if (table_array == NULL) return NULL;
  Py_DECREF(table_obj);
  table = (double*)table_array->data;
  kind = table[0];
  version = table[1];
  if (kind != 5 || version != 2) {
    PyErr_Format (PyExc_ValueError, "Only support kind = 5, version = 2.  Found: kind = %i, version = %i", kind, version);
    return NULL;
  }
  ptop = table[3];
  pref = table[4];
  rcoef1 = table[5];
  rcoef2 = table[6];
  {
    unsigned long long x = *(unsigned long long*)(table+7);
    ref_name[0] = x & 255;
    ref_name[1] = (x>>8) & 255;
    ref_name[2] = (x>>16) & 255;
    ref_name[3] = x>>24;
  }

  // Get number of thermodynamic & momentum levels.
  double *tab_m, *tab_t;
  {
    int skip = table[2];
    int nj = PyArray_DIMS(table_array)[1];
    int nk = (nj-3-skip)/2;
    n_m = nk+1;
    n_t = nk+2;
    tab_m = table + 3*skip;
    tab_t = table + 3*(skip+nk+1);
  }
  // Allocate output arrays
  IP1_m = (PyArrayObject*) PyArray_SimpleNew(1, &n_m, NPY_INT);
  A_m = (PyArrayObject*) PyArray_SimpleNew(1, &n_m, NPY_FLOAT64);
  B_m = (PyArrayObject*) PyArray_SimpleNew(1, &n_m, NPY_FLOAT64);
  IP1_t = (PyArrayObject*) PyArray_SimpleNew(1, &n_t, NPY_INT);
  A_t = (PyArrayObject*) PyArray_SimpleNew(1, &n_t, NPY_FLOAT64);
  B_t = (PyArrayObject*) PyArray_SimpleNew(1, &n_t, NPY_FLOAT64);
  ip1_m = (int*)PyArray_DATA(IP1_m);
  a_m = (double*)PyArray_DATA(A_m);
  b_m = (double*)PyArray_DATA(B_m);
  ip1_t = (int*)PyArray_DATA(IP1_t);
  a_t = (double*)PyArray_DATA(A_t);
  b_t = (double*)PyArray_DATA(B_t);

  // Read vertical information from the table
  for (i = 0; i < n_m; i++) {
    ip1_m[i] = tab_m[0];
    a_m[i] = tab_m[1];
    b_m[i] = tab_m[2];
    tab_m+=3;
  }
  for (i = 0; i < n_t; i++) {
    ip1_t[i] = tab_t[0];
    a_t[i] = tab_t[1];
    b_t[i] = tab_t[2];
    tab_t+=3;
  }
  Py_DECREF(table_array);

  PyObject *ret = Py_BuildValue ("iiffffs#OOOOOO", kind, version, ptop, pref, rcoef1, rcoef2, ref_name, 4, IP1_m, A_m, B_m, IP1_t, A_t, B_t);
  Py_DECREF(IP1_m);
  Py_DECREF(A_m);
  Py_DECREF(B_m);
  Py_DECREF(IP1_t);
  Py_DECREF(A_t);
  Py_DECREF(B_t);

  return ret;
}

// Find specific A and B from the given table
static PyObject *get_loghybrid_a_b (PyObject *self, PyObject *args) {
  PyArrayObject *IP1, *A, *B, *IP1_m, *A_m, *B_m, *IP1_t, *A_t, *B_t;
  int *ip1, *ip1_m, *ip1_t;
  double *a, *b, *a_m, *b_m, *a_t, *b_t;
  npy_intp n, n_m, n_t;
  int i, j;

  // Parse inputs
  if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!", &PyArray_Type, &IP1, &PyArray_Type, &IP1_m, &PyArray_Type, &A_m, &PyArray_Type, &B_m, &PyArray_Type, &IP1_t, &PyArray_Type, &A_t, &PyArray_Type, &B_t)) return NULL;
  // Assume arrays are of the right type and contiguous (and of the right length).
  n_m = PyArray_SIZE(IP1_m);
  n_t = PyArray_SIZE(IP1_t);
  n = PyArray_SIZE(IP1);

  // Define output arrays
  A = (PyArrayObject*)PyArray_SimpleNew(1, &n, NPY_FLOAT64);
  B = (PyArrayObject*)PyArray_SimpleNew(1, &n, NPY_FLOAT64);

  // Match the levels
  ip1 = (int*)PyArray_DATA(IP1);
  a = (double*)PyArray_DATA(A);
  b = (double*)PyArray_DATA(B);
  ip1_m = (int*)PyArray_DATA(IP1_m);
  a_m = (double*)PyArray_DATA(A_m);
  b_m = (double*)PyArray_DATA(B_m);
  ip1_t = (int*)PyArray_DATA(IP1_t);
  a_t = (double*)PyArray_DATA(A_t);
  b_t = (double*)PyArray_DATA(B_t);
  for (i = 0; i < n; i++) {
    a[i] = b[i] = 0.;
    for (j = 0; j < n_m; j++) {
      if (ip1[i] == ip1_m[j]) {
        a[i] = a_m[j];
        b[i] = b_m[j];
        break;
      }
    }
    if (j < n_m) continue;
    for (j = 0; j < n_t; j++) {
      if (ip1[i] == ip1_t[j]) {
        a[i] = a_t[j];
        b[i] = b_t[j];
        break;
      }
    }
  }

  PyObject *ret = Py_BuildValue ("OO", A, B);
  Py_DECREF(A);
  Py_DECREF(B);
  return ret;
}


// make_hy_record(eta.values,eta.A,eta.B)
static PyObject *make_hy_record (PyObject *self, PyObject *args) {
  float ptop, rcoef, pref;

  if (!PyArg_ParseTuple(args, "fff", &ptop, &rcoef, &pref)) return NULL;

  // Encode in a record
  Py_INCREF(descr);
  npy_intp dims[] = {1};
  PyArrayObject *record = (PyArrayObject*) PyArray_SimpleNewFromDescr (1, dims, descr);
  if (record == NULL) return NULL;

  HEADER *hy = (HEADER*)record->data;

  // Encode IG1, IG2, IG3, IG4
  {
    float x1, x2;
    f77name(hybref_to_ig)(&hy->ig1, &hy->ig2, &hy->ig3, &hy->ig4, &rcoef, &pref, &x1, &x2);
  }

  // Encode ptop
  {
    int kind = 2, mode = 2, flag = 0;
    f77name(convip)(&hy->ip1, &ptop, &kind, &mode, "", &flag);
  }

  // Other parameters
  strncpy (hy->nomvar, "HY  ", 4);
  strncpy (hy->grtyp, "X ", 2);
  strncpy (hy->typvar, "X ", 2);
  hy->datyp = 5;
  hy->nbits = 32;
  hy->ni = 1;
  hy->nj = 1;
  hy->nk = 1;

  // NOTE: data function needs to be set from the Python side

  return (PyObject*)record;

}

// Helper methods - find a coordinate record for the given field
// (match ig1,ig2,ig3 to a coordinate's ip1,ip2,ip3).
int find_coord(HEADER *records, int num_records, char nomvar[4], int ig1, int ig2, int ig3, int ni, int nj) {
  int i;
  for (i = 0; i < num_records; i++) {
    if (strncmp(records[i].nomvar,nomvar,4) != 0) continue;
    if (records[i].ni != ni) continue;
    if (records[i].nj != nj) continue;
    if (records[i].ip1 != ig1) continue;
    if (records[i].ip2 != ig2) continue;
    if (records[i].ip3 != ig3) continue;
    return i;
  }
  return -1;
}
int find_xrec(HEADER *records, int n, int varid) {
  char grtyp = records[varid].grtyp[0];
  int ig1 = records[varid].ig1;
  int ig2 = records[varid].ig2;
  int ig3 = records[varid].ig3;
  int ni = records[varid].ni;
  int nj = records[varid].nj;
  if (grtyp == 'Y')
    return find_coord(records, n, ">>  ", ig1, ig2, ig3, ni, nj);
  if (grtyp == 'Z')
    return find_coord(records, n, ">>  ", ig1, ig2, ig3, ni, 1);
  return -1;
}
int find_yrec(HEADER *records, int n, int varid) {
  char grtyp = records[varid].grtyp[0];
  int ig1 = records[varid].ig1;
  int ig2 = records[varid].ig2;
  int ig3 = records[varid].ig3;
  int ni = records[varid].ni;
  int nj = records[varid].nj;
  if (grtyp == 'Y')
    return find_coord(records, n, "^^  ", ig1, ig2, ig3, ni, nj);
  if (grtyp == 'Z')
    return find_coord(records, n, "^^  ", ig1, ig2, ig3, 1, nj);
  return -1;
}


// Helper functions - determine if latitude/longitude fields are really 1D
int lat_is_1d (float *lat, int ni, int nj) {
  int i, j;
  for (j = 0; j < nj; j++, lat+=ni) {
    for (i = 1; i < ni; i++) {
      if (fabsf(lat[0]-lat[i]) > 1E-4) return 0;
    }
  }
  return 1;
}
int lon_is_1d (float *lon, int ni, int nj) {
  int i, j;
  float *ref = lon;
  for (j = 0; j < nj; j++, lon+=ni) {
    for (i = 0; i < ni; i++) {
      if (fabsf(ref[i]-lon[i]) > 1E-4) return 0;
    }
  }
  return 1;
}

// Construct latitude/longitude fields from the given records
static PyObject *get_latlon (PyObject *self, PyObject *args) {

  // Get the input records
  PyArrayObject *record_array;
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &record_array)) return NULL;
  if (record_array->descr != descr) return NULL;
  if (PyArray_ITEMSIZE(record_array) != sizeof(HEADER)) return NULL;
  record_array = PyArray_GETCONTIGUOUS(record_array);
  HEADER *records = (HEADER*)record_array->data;
  int num_records = PyArray_SIZE(record_array);

  // Create a dictionary to hold unique latitude & longitude arrays.
  // Keys are (grtyp,ig1,ig2,ig3,ig4,ni,nj)
  PyObject *dict = PyDict_New();

  int i;
  for (i = 0; i < num_records; i++) {
    // Skip coordinate records
    if (strncmp(records[i].nomvar,">>  ",4) == 0) continue;
    if (strncmp(records[i].nomvar,"^^  ",4) == 0) continue;
    if (strncmp(records[i].nomvar,"HY  ",4) == 0) continue;
    if (strncmp(records[i].nomvar,"!!  ",4) == 0) continue;
    // Construct the key for this record
    char *grtyp = records[i].grtyp;
    int ig1 = records[i].ig1;
    int ig2 = records[i].ig2;
    int ig3 = records[i].ig3;
    int ig4 = records[i].ig4;
    int ni = records[i].ni;
    int nj = records[i].nj;
    PyObject *key = Py_BuildValue("s#iiiiii",grtyp,1,ig1,ig2,ig3,ig4,ni,nj);
    if (key == NULL) return NULL;
    // Check if we've already handled this grid
    if (PyDict_Contains(dict,key)) continue;

    // Needed for error messages
    char nomvar[] = "    ";
    strncpy (nomvar, records[i].nomvar, 4);

    // Get the grid id
    int gdid;
    switch (*grtyp) {
      case 'Z':;
        int xrec = find_xrec (records, num_records, i);
        int yrec = find_yrec (records, num_records, i);

        if (xrec < 0 || yrec < 0) {
          PyErr_Format (PyExc_ValueError, "Coordinate record(s) not found for '%s'", nomvar);
          return NULL;
        }

        PyArrayObject *xobj = (PyArrayObject*)PyObject_CallObject(records[xrec].data_func,NULL);
        PyArrayObject *yobj = (PyArrayObject*)PyObject_CallObject(records[yrec].data_func,NULL);
        float *ax = (float*)xobj->data;
        float *ay = (float*)yobj->data;
        char *grref = records[xrec].grtyp;
        int ig1ref = records[xrec].ig1;
        int ig2ref = records[xrec].ig2;
        int ig3ref = records[xrec].ig3;
        int ig4ref = records[xrec].ig4;
        gdid = c_ezgdef_fmem(ni, nj, grtyp, grref, 
                ig1ref, ig2ref, ig3ref, ig4ref, ax, ay);
        Py_DECREF(xobj);
        Py_DECREF(yobj);
        break;
      case 'A':
      case 'B':
      case 'E':
      case 'G':
      case 'L':
      case 'N':
      case 'S':
        gdid = c_ezqkdef (ni, nj, grtyp, ig1, ig2, ig3, ig4, 0);
        break;
      default:
        PyErr_Format (PyExc_ValueError, "Unhandled grid type '%c'", *grtyp);
        return NULL;
    }

    // Extract the latitudes and longitudes
    npy_intp dims[] = {nj,ni};
    PyArrayObject *lat = (PyArrayObject*)PyArray_SimpleNew (2, dims, NPY_FLOAT32);
    PyArrayObject *lon = (PyArrayObject*)PyArray_SimpleNew (2, dims, NPY_FLOAT32);
    c_gdll (gdid, (float*)lat->data, (float*)lon->data);

    // Extract x and y coordinates
    PyArrayObject *ax = (PyArrayObject*)PyArray_SimpleNew (1, dims+1, NPY_FLOAT32);
    PyArrayObject *ay = (PyArrayObject*)PyArray_SimpleNew (1, dims+0, NPY_FLOAT32);
    c_gdgaxes (gdid, (float*)ax->data, (float*)ay->data);

    // Done with the grid
    c_gdrls (gdid);


    // Can the latitudes and longitudes be reduced to 1D arrays?
    if (lon_is_1d ((float*)lon->data, ni, nj)) {
//      printf ("1D longitudes detected\n");
      PyArrayObject *newlon = (PyArrayObject*)PyArray_SimpleNew (1, dims+1, NPY_FLOAT32);
      memcpy (newlon->data, lon->data, ni*sizeof(float));
      Py_DECREF(lon);
      lon = newlon;
    }
    if (lat_is_1d ((float*)lat->data, ni, nj)) {
//      printf ("1D latitudes detected\n");
      PyArrayObject *newlat = (PyArrayObject*)PyArray_SimpleNew (1, dims+0, NPY_FLOAT32);
      int j;
      for (j = 0; j < nj; j++) {
        ((float*)newlat->data)[j] = ((float*)lat->data)[j*ni];
      }
      Py_DECREF(lat);
      lat = newlat;
    }

    // Build the value
    PyObject *value = Py_BuildValue("OOOO", ax, ay, lat, lon);
    Py_DECREF(ax);
    Py_DECREF(ay);
    Py_DECREF(lat);
    Py_DECREF(lon);
    PyDict_SetItem(dict,key,value);
    Py_DECREF(key);
    Py_DECREF(value);
  }

  // Clean up local references
  Py_DECREF(record_array);

  // Return the arrays
  return dict;
}

static PyMethodDef FST_Methods[] = {
  {"read_records", fstd_read_records, METH_VARARGS, "Get all record headers from an FSTD file"},
  {"write_records", fstd_write_records, METH_VARARGS, "Write a set of records into a given FSTD file"},
  {"stamp2date", stamp2date, METH_VARARGS, "Convert CMC timestamps to seconds since 1980-01-01 00:00:00"},
  {"date2stamp", date2stamp, METH_VARARGS, "Convert seconds since 1980-01-01 00:00:00 to a CMC timestamp"},
  {"decode_levels", decode_levels, METH_VARARGS, "Decode vertical levels"},
  {"encode_levels", encode_levels, METH_VARARGS, "Encode vertical levels"},
  {"get_hybrid_a_b", get_hybrid_a_b, METH_VARARGS, "Get A and B arrays from HY record and specified levels"},
  {"get_loghybrid_table", get_loghybrid_table, METH_VARARGS, "Get info table from !! record"},
  {"get_loghybrid_a_b", get_loghybrid_a_b, METH_VARARGS, "Get A and B from table and specific ip1 values"},
  {"make_hy_record", make_hy_record, METH_VARARGS, "Construct an HY record from the given eta, a, and b arrays"},
  {"get_latlon", get_latlon, METH_VARARGS, "Create lat/lon arrays from the given records"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initfstd_core(void) {
  PyObject *m = Py_InitModule("fstd_core", FST_Methods);
  import_array();

  Py_INCREF(&FSTD_Unit_Type);
  if (PyType_Ready(&FSTD_Unit_Type) < 0) return;
  PyModule_AddObject (m, "FSTD_Unit", (PyObject*)&FSTD_Unit_Type);

  Py_INCREF(&RecordGetter_Type);
  if (PyType_Ready(&RecordGetter_Type) < 0) return;
  PyModule_AddObject (m, "RecordGetter", (PyObject*)&RecordGetter_Type);

  PyObject *header_structure = Py_BuildValue("[(s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s), (s,s)]", "pad", "i4", "dateo", "i4", "deet", "i4", "npas", "i4", "ni", "i4", "nj", "i4", "nk", "i4", "nbits", "i4", "datyp", "i4", "ip1", "i4", "ip2", "i4", "ip3", "i4", "typvar", "a2", "nomvar", "a4", "etiket", "a12", "grtyp", "a2", "ig1", "i4", "ig2", "i4", "ig3", "i4", "ig4", "i4", "swa", "i4", "lng", "i4", "dltf", "i4", "ubc", "i4", "extra1", "i4", "extra2", "i4", "extra3", "i4", "data_func", "O");
  if (header_structure == NULL) return;
  PyArray_DescrConverter (header_structure, &descr);
  Py_DECREF(header_structure);
  PyModule_AddObject (m, "record_descr", (PyObject*)descr);

}

