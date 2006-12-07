/*******************************************
** macro_bunch_store.h
** Contains:
** 
*******************************************/

#ifndef HAVE_MACRO_BUNCH_STORE_H
#define HAVE_MACRO_BUNCH_STORE_H

#include <iostream>

#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>

using namespace boost::python;

struct Macro_bunch_store
{
  numeric::array* numeric_local_particles;
  double* local_particles;
  int local_num, total_num;
  double total_current;
  numeric::array* numeric_units;
  double* units;
  numeric::array* numeric_ref_particle;
  double* ref_particle;
  bool is_z;

  Macro_bunch_store(numeric::array& numeric_local_particles,
	      int local_num, int total_num, double total_current,
	      numeric::array& numeric_units,
	      numeric::array& numeric_ref_particle,
	      bool is_z);
  ~Macro_bunch_store();
};

#endif				//	ifndef  HAVE_SCALAR_FIELD_H


