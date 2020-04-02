#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>

#include "portable_import_array.h"
#include "numpy/ndarraytypes.h"
#include "numpy/__multiarray_api.h"

void portable_import_array() {
  if (_import_array() != 0) {
    PyErr_Print();
    PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
  }
  return;
}

