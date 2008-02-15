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
#include <string>

#include "array_nd/array_nd.h"
#include "array_nd/array_2d.h"
#include "array_nd/array_1d.h"
#include "array_nd/array_nd_python.h"
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

void
chef_propagate(numeric::array& numeric_particles, int num_particles,
	       bmlnElmnt& element, double energy_in, std::string particle_type,
	       numeric::array& numeric_u_in, numeric::array& numeric_u_out)
{
    Array_1d<double> u_in =
        Array_nd_from_PyObject<double>(numeric_u_in.ptr());
    Array_1d<double> u_out =
        Array_nd_from_PyObject<double>(numeric_u_out.ptr());
        
    Array_2d<double> particles = 
        Array_nd_from_PyObject<double>(numeric_particles.ptr());

  Vector chef_state(6);
  for(int part=0; part<num_particles; ++part) {
    for (int impact_index=0; impact_index<6; ++impact_index) {
      int chef_index = convert_chef_index(impact_index);
      chef_state[chef_index] = particles(impact_index,part)/
	u_in(impact_index);
    }
    Particle *particle_ptr;
    if (particle_type == "proton") {
        particle_ptr = new Proton(energy_in);
    } else if (particle_type == "positron") {
        particle_ptr = new Positron(energy_in);
    } else {
        throw 
            std::runtime_error("chef_propagate: unknown particle_type " + particle_type);
    }
    particle_ptr->State() = chef_state;
    element.propagate(*particle_ptr);
    chef_state = particle_ptr->State();
    for (int impact_index=0; impact_index<6; ++impact_index) {
      int chef_index = convert_chef_index(impact_index);
      particles(impact_index,part) = chef_state[chef_index]*
	u_out(impact_index);
    }
    delete particle_ptr;
  }
}

    
    
BOOST_PYTHON_MODULE(chef_propagate)
{
  numeric::array::set_module_and_type("Numeric", "ArrayType");
  def("chef_propagate",&chef_propagate);
}

