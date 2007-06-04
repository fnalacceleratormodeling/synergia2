#include <iostream>
#include <list>
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <Numeric/arrayobject.h>
#include <vector>
#include "bmlfactory/bmlfactory.h"
#include "mxyzptlk/Mapping.h"
#include "beamline/bmlnElmnt.h"
#include "beamline/Particle.h"

extern "C" {
#include <sys/time.h>
}

double 
double_time()
{
  timeval t;
  gettimeofday(&t,NULL);
  return t.tv_sec + t.tv_usec/1.0e6;
}
  
using namespace boost::python;

int convert_chef_index(int impact_index)
{
  return impact_index/2+3*(impact_index%2);
}

class Particles
{
private:
  double *data;
  int index(int row,int col) {
    return col*7 + row;
  };
public:
  Particles(numeric::array& numeric_particles) {
    data = reinterpret_cast<double *>
      (reinterpret_cast<PyArrayObject*>(numeric_particles.ptr())->data);
  };
  double & operator()(int component, int particle) {
    return data[index(component,particle)];
  };
};


class Unit_conversion
{
private:
  double *data;
public:
  Unit_conversion(numeric::array& numeric_u) {
    data = reinterpret_cast<double *>
    (reinterpret_cast<PyArrayObject*>(numeric_u.ptr())->data);
  };
  double & operator()(int impact_index) {
    return data[impact_index];
  };
};

void
chef_propagate(numeric::array& numeric_particles, int num_particles,
	       bmlnElmnt& element, double energy_in,
	       numeric::array& numeric_u_in, numeric::array& numeric_u_out)
{
  Unit_conversion u_in(numeric_u_in);
  Unit_conversion u_out(numeric_u_out);

  Particles particles(numeric_particles);

  Vector chef_state(6);
  for(int part=0; part<num_particles; ++part) {
    for (int impact_index=0; impact_index<6; ++impact_index) {
      int chef_index = convert_chef_index(impact_index);
      chef_state[chef_index] = particles(impact_index,part)/
	u_in(impact_index);
    }
    Proton proton(energy_in);
    proton.setState(chef_state);
    element.propagate(proton);
    chef_state = proton.getState();
    for (int impact_index=0; impact_index<6; ++impact_index) {
      int chef_index = convert_chef_index(impact_index);
      particles(impact_index,part) = chef_state[chef_index]*
	u_out(impact_index);
    }
  }
}

    
    
BOOST_PYTHON_MODULE(chef_propagate)
{
  numeric::array::set_module_and_type("Numeric", "ArrayType");
  def("chef_propagate",&chef_propagate);
}

