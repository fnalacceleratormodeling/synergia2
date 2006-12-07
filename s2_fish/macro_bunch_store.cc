#include "macro_bunch_store.h"
#include <Numeric/arrayobject.h>

using namespace boost::python;

template<class T> T data_from_numeric_array(numeric::array& array)
{
  return reinterpret_cast<T>
    (reinterpret_cast<PyArrayObject*>(array.ptr())->data);
}

Macro_bunch_store::Macro_bunch_store(numeric::array& numeric_local_particles,
			 int local_num, int total_num, double total_current,
			 numeric::array& numeric_units,
			 numeric::array& numeric_ref_particle,
			 bool is_z)
{
  this->numeric_local_particles = &numeric_local_particles;
  local_particles = data_from_numeric_array<double *>(numeric_local_particles);
  this->local_num = local_num;
  this->total_num = total_num;
  this->total_current = total_current;
  this->numeric_units = &numeric_units;
  units = data_from_numeric_array<double*>(numeric_units);
  this->numeric_ref_particle = &numeric_ref_particle;
  ref_particle = data_from_numeric_array<double*>(numeric_ref_particle);
  this->is_z = is_z;
}

Macro_bunch_store::~Macro_bunch_store()
{
}
