#include "ltwt_containers.h"
#include <Numeric/arrayobject.h>

Ltwt_array::Ltwt_array(numeric::array& numeric_array)
{
  this->numeric_array = &numeric_array;
  PyArrayObject* array_object = reinterpret_cast<PyArrayObject*>
    (numeric_array.ptr());
  data = reinterpret_cast<double*>(array_object->data);
  stride = array_object->strides[0]/sizeof(double);
}

double& Ltwt_array::operator()(int index)
{
  return *(data + stride*index);
}

Ltwt_matrix::Ltwt_matrix(numeric::array& numeric_array)
{
  this->numeric_array = &numeric_array;
  PyArrayObject* array_object = reinterpret_cast<PyArrayObject*>
    (numeric_array.ptr());
  data = reinterpret_cast<double*>(array_object->data);
  stride0 = array_object->strides[0]/sizeof(double);
  stride1 = array_object->strides[1]/sizeof(double);
}

double& Ltwt_matrix::operator()(int row, int col)
{
  return *(data + stride0*row + stride1*col);
}

  
