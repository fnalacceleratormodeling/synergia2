#include <iostream>
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <Numeric/arrayobject.h>
#include <vector>
#include "bmlfactory/bmlfactory.h"

class Term
{
private:
  int order;
public:
  Term(int order_in) { order = order_in; };
};

using namespace boost::python;

// indexing for particles
int pindex(int row,int col)
{
  return col*7 + row;
}

// indexing for maps
int mindex(int row,int col)
{
  return row*7 + col;
}

void apply_map1(numeric::array& numeric_particles, int num_particles,
		numeric::array& numeric_map)
{
  double *particles;
  particles = reinterpret_cast<double *>
    (reinterpret_cast<PyArrayObject*>(numeric_particles.ptr())->data);
  double *map;
  map = reinterpret_cast<double *>
    (reinterpret_cast<PyArrayObject*>(numeric_map.ptr())->data);
  
  double temp[6];
  for(int part=0; part<num_particles; ++part) {
    for(int i=0; i<6; ++i) {
      temp[i] = 0.0;
      for(int j=0; j<6; ++j) {
	temp[i] += particles[pindex(j,part)]*map[mindex(i,j)];
      }
    }
    for(int i=0; i<6; ++i) {
      particles[pindex(i,part)] = temp[i];
    }
  }
}
	
BOOST_PYTHON_MODULE(apply_map)
{
  numeric::array::set_module_and_type("Numeric", "ArrayType");
  def("apply_map1",&apply_map1);
}

