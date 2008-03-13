// -*- C++ -*-
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

#include "array_nd/array_2d.h"
#include "array_nd/array_1d.h"

using namespace boost::python;

struct Macro_bunch_store {
    Array_2d<double> local_particles;
    int local_num, total_num;
    double mass;
    int charge;
    double total_current;
    Array_1d<double> units; // Unit conversion: X^impact_i = units_i X^real_i
    Array_1d<double> ref_particle;
    numeric::array *numeric_local_particles, *numeric_units,
    *numeric_ref_particle;
    bool is_fixedz;

    Macro_bunch_store(numeric::array& numeric_local_particles,
                      int local_num, int total_num,
                      double mass, int charge, double total_current,
                      numeric::array& numeric_units,
                      numeric::array& numeric_ref_particle,
                      bool is_fixedz);
    numeric::array get_local_particles();
    numeric::array get_units();
    numeric::array get_ref_particle();
    void convert_to_fixedt();
    void convert_to_fixedz();
    void check(numeric::array & array);
    double get_coord(int coord_index, int particle_index);
    void set_local_particles(numeric::array &particles);
    ~Macro_bunch_store();
};

#endif				//	ifndef  HAVE_MACRO_BUNCH_STORE_H


