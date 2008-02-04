#define NUMARRAY youbetcha

#include <iostream>
#include <octave/config.h>
#include <octave/octave.h>
#include "symtab.h"
#include "parse.h"
#include "unwind-prot.h"
#include "toplev.h"
#include "error.h" 
#include "quit.h"
#include "variables.h"
#include "sighandlers.h"
#include "sysdep.h"
#include <octave/oct.h>
#include <octave/symtab.h>
#include <string>

#include "octave_class.h"

Octave::Octave()
{
#ifdef NUMARRAY
import_libnumarray();
#endif
  import_array();
  int argc = 2;
  char * argv [] = {"octapy","-q"};
  octave_main(argc,argv,1);
}

Octave::~Octave()
{
  do_octave_atexit();
}

static void
recover_from_exception (void)
{
  unwind_protect::run_all ();
  can_interrupt = true;
  octave_interrupt_immediately = 0;
  octave_interrupt_state = 0;
  octave_allocation_error = 0;
  octave_restore_signal_mask ();
  octave_catch_interrupts ();
}

int 
Octave::execute(std::string string)
{
  int parse_status;

  octave_save_signal_mask ();
  if (octave_set_current_context)
    {
#if defined (USE_EXCEPTIONS_FOR_INTERRUPTS)
      panic_impossible ();
#else
      unwind_protect::run_all ();
      raw_mode (0);
      std::cout << "\n";
      octave_restore_signal_mask ();
#endif
    }

  can_interrupt = true;
  octave_catch_interrupts ();
  octave_initialized = true;

  try 
    {
      curr_sym_tab = top_level_sym_tab;
      reset_error_handler ();
      eval_string(string, false, parse_status);
    }
  catch (octave_interrupt_exception)
    {
      recover_from_exception ();
      std::cout << "\n"; 
      error_state = -2; 
    }
  catch (std::bad_alloc)
    {
      recover_from_exception ();
      std::cout << "\n"; 
      error_state = -3;
    }

  octave_restore_signal_mask();
  octave_initialized = false;

  // XXX FIXME XXX callbacks calling embed_octave
  // may or may not want error_state reset.
  return error_state;
}

PyObject *
Octave::top_level_vars()
{
  PyObject * result = PyList_New(0);
  string_vector name_list = top_level_sym_tab->variable_name_list();
  for(int i=0; i <name_list.length() ; i++) {
    PyList_Append(result,PyString_FromString(name_list[i].c_str()));
  }
  name_list = top_level_sym_tab->global_variable_name_list();
  for(int i=0; i <name_list.length() ; i++) {
    PyList_Append(result,PyString_FromString(name_list[i].c_str()));
  }
  return result;
}

octave_value
get_top_level_value (const std::string& nm)
{
  octave_value retval;
  symbol_record *sr = top_level_sym_tab->lookup (nm);
  if (sr) {
    octave_value val = sr->variable_value ();
    if (val.is_undefined ())
      error ("get_top_level_value: undefined symbol `%s'", nm.c_str ());
    else
      retval = val;
  } else {
    error ("get_top_level_value: unknown symbol `%s'", nm.c_str ());
  }
  return retval;
}

void
set_top_level_value (const std::string& nm, octave_value & value)
{
  symbol_record *sr = top_level_sym_tab->lookup (nm, true);
  if (sr){
    sr->define(value);
  }
  else
    error ("set_top_level_value: unknown symbol `%s'", nm.c_str ());
}

PyObject *
Octave::get_value(std::string name)
{
  PyObject *value = NULL;
  octave_value oct_val = get_top_level_value(name);
  if (oct_val.is_scalar_type()) {
    if (oct_val.is_real_type()) {
      value = Py_BuildValue("d",oct_val.double_value());
    } else {
      value = PyComplex_FromDoubles(oct_val.complex_value().real(),
				    oct_val.complex_value().imag());
    }
  } else {
    if (oct_val.is_string()) {
      value = Py_BuildValue("s",oct_val.string_value().c_str());
    } else {
      if (oct_val.is_real_matrix()) {
	int rank = oct_val.matrix_value().dims().length();
	int dims[rank], vector_dims[1], matrix_dims[2];
	bool is_vector = false;
	int vector_length = 1;
	for (int i = 0; i< rank ; ++i) {
	  dims[i] = oct_val.matrix_value().dims()(i);
	  if (rank == 2 && dims[i] == 1) {
	    is_vector = true;
	  } else {
	    vector_length = dims[i];
	  }
	}
	if (is_vector) {
	  vector_dims[0] = vector_length;
#ifdef NUMARRAY
	  value = (PyObject *) NA_vNewArray(oct_val.matrix_value().fortran_vec(),
					    tFloat64,1,vector_dims);
#else
	  value = PyArray_FromDimsAndData(1,vector_dims,PyArray_DOUBLE,
					  (char*)(oct_val.matrix_value().fortran_vec()));
#endif
	} else {
	  matrix_dims[0] = dims[1];
	  matrix_dims[1] = dims[0];
#ifdef NUMARRAY
	  value = (PyObject *) 
	    NA_vNewArray(oct_val.matrix_value().fortran_vec(),
			 tFloat64,rank,matrix_dims);
#else
	  value = 
	    PyArray_FromDimsAndData(2,matrix_dims,PyArray_DOUBLE,
				    (char*)(oct_val.matrix_value().fortran_vec()));
#endif
	}
      }
    }
  }
  if (!value) {
    std::cout << "Octave.get_value: " << name << " has type " 
	      << oct_val.type_name() << ", which is not handled" << std::endl;
    value = Py_BuildValue("");
  }
  return value;
}

void
Octave::set_value(std::string name, PyObject * value)
{
  octave_value * oct_val = NULL;
  if (PyInt_Check(value)) {
    long int cxx_value(PyInt_AsLong(value));
    oct_val = new octave_value(cxx_value);
  } else if (PyFloat_Check(value)) {
    double cxx_value(PyFloat_AsDouble(value));
    oct_val = new octave_value(cxx_value);
  } else if (PyComplex_Check(value)) {
    Complex cxx_value(PyComplex_RealAsDouble(value),
			   PyComplex_ImagAsDouble(value));
    oct_val = new octave_value(cxx_value);
  } else if (PyString_Check(value)) {
    std::string cxx_value(PyString_AsString(value));
    oct_val = new octave_value(cxx_value);
  } else if (PyArray_Check(value)) {
#ifdef NUMARRAY
    PyArrayObject * pa = NA_InputArray(value,tFloat64,NUM_C_ARRAY);
#else
    PyArrayObject * pa = (PyArrayObject *)value;
#endif
    if (pa->nd == 1) {
      ColumnVector cv(pa->dimensions[0]);
      for (int i=0; i<pa->dimensions[0]; ++i) {
	cv(i) = ((double *) pa->data)[i];
      }
      oct_val = new octave_value(cv);
    } else if (pa->nd == 2) {
      Matrix m(pa->dimensions[0],pa->dimensions[1]);
      for (int i=0; i<pa->dimensions[0]; ++i) {
	for (int j=0; j<pa->dimensions[1]; ++j) {
	  m(i,j) = ((double *) pa->data)[i*pa->dimensions[1] + j];
	}
      }
      oct_val = new octave_value(m);
    } else {
      std::cout << "Arrays of dimension " << pa->nd << " are not supported.\n";
    }
  }
			   
  if (oct_val) {
    set_top_level_value(name,(*oct_val));
  } else {
    std::cout << "Octave:set_value: variable type not handled.\n";
  }
}
