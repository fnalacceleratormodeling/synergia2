/*******************************************
** ltwt_containers.h
** Contains:
** 
*******************************************/

#ifndef HAVE_LTWT_CONTAINERS_H
#define HAVE_LTWT_CONTAINERS_H

#include <iostream>

#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>

using namespace boost::python;

struct Ltwt_array
{
  numeric::array* numeric_array;
  double* data;
  size_t stride;
  Ltwt_array(numeric::array& numeric_array);
  double& operator()(int index);
};

struct Ltwt_matrix
{
  numeric::array* numeric_array;
  double* data;
  size_t stride0,stride1;
  Ltwt_matrix(numeric::array& numeric_array);
  double& operator()(int row, int col);
};

#endif // HAVE_LTWT_CONTAINERS_H
