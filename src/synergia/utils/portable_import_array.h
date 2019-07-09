#ifndef PORTABLE_IMPORT_ARRAY_H_
#define PORTABLE_IMPORT_ARRAY_H_

// This function is here to deal with converting from Python2 to Python3.
// In Python2, the import_array() macro returns nothing; in Python3, it
// returns NULL.
// The BOOST_PYTHON_MODULE() macro is using it in a function that returns
// nothing, so the Python3 implementation causes compilation failures.
//
// This function (not macro) reproduces the implementation of the
// import_array() macro.
//
inline
void portable_import_array() {
  if (_import_array() < 0) {
    PyErr_Print();
    PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
  }
  return;
}

#endif /* PORTABLE_IMPORT_ARRAY_H_ */
