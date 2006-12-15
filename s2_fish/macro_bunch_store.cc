#include "macro_bunch_store.h"
#include <Numeric/arrayobject.h>

using namespace boost::python;

template<class T> T data_from_numeric_array(numeric::array& array)
{
  return reinterpret_cast<T>
    (reinterpret_cast<PyArrayObject*>(array.ptr())->data);
}

Macro_bunch_store::Macro_bunch_store(numeric::array& numeric_local_particles,
				     int local_num, int total_num, 
				     double total_current,
				     numeric::array& numeric_units,
				     numeric::array& numeric_ref_particle,
				     bool is_fixedz) :
  local_particles(numeric_local_particles),
  units(numeric_units),
  ref_particle(numeric_ref_particle)
{
  this->local_num = local_num;
  this->total_num = total_num;
  this->total_current = total_current;
  this->is_fixedz = is_fixedz;
}

numeric::array Macro_bunch_store::get_local_particles()
{
  return *(local_particles.numeric_array);
}

numeric::array Macro_bunch_store::get_units()
{
  return *(units.numeric_array);
}

numeric::array Macro_bunch_store::get_ref_particle()
{
  return *(ref_particle.numeric_array);
}

void Macro_bunch_store::convert_to_fixedt()
{
  if (is_fixedz) {
    double gamma = -ref_particle(5);
    for(int i=0; i<local_num; ++i) {
      local_particles(0,i) /= units(0);
      local_particles(2,i) /= units(2);
      double xp = local_particles(1,i);
      double yp = local_particles(3,i);
      double rcp_gammai = 1.0/(gamma - local_particles(5,i));
      double betai = sqrt(1.0 - rcp_gammai*rcp_gammai * (1 + xp*xp +yp*yp));
      local_particles(4,i) *= -gamma*betai/units(0); // units(0)
                                                     // is not
                                                     // an error!
    }
    is_fixedz = false;
  }
}

void Macro_bunch_store::convert_to_fixedz()
{
  if (! is_fixedz) {
    double gamma = -ref_particle(5);
    for(int i=0; i<local_num; ++i) {
      local_particles(0,i) *= units(0);
      local_particles(2,i) *= units(2);
      double xp = local_particles(1,i);
      double yp = local_particles(3,i);
      double rcp_gammai = 1.0/(gamma - local_particles(5,i));
      double betai = sqrt(1.0 - rcp_gammai*rcp_gammai * (1 + xp*xp +yp*yp));
      local_particles(4,i) /= -gamma*betai/units(0); // units(0)
					  	     // is not
						     // an error!
    }
    is_fixedz = true;
  }
}

double Macro_bunch_store::get_coord(int coord_index,int particle_index)
{
  return local_particles(coord_index,particle_index);
}

Macro_bunch_store::~Macro_bunch_store()
{
}
