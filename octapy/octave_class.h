#ifndef OCTAVE_CLASS_H__
#define OCTAVE_CLASS_H__
#undef HAVE_STAT
#undef _POSIX_C_SOURCE
#include <Python.h>
#include <string>
#ifdef NUMARRAY
#include <numarray/libnumarray.h>
#include <numarray/arrayobject.h>
#else
#include <numpy/arrayobject.h>
#endif
class Octave
{
 public:
  Octave();
  int execute(std::string string);
  PyObject * top_level_vars();
  PyObject * get_value(std::string name);
  void set_value(std::string name, PyObject *value);
  ~Octave();
};


#endif

