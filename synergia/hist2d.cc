#include <iostream>
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <numpy/arrayobject.h>

using namespace boost::python;

void
hist2d(numeric::array& numeric_x, double xmin, double xmax, long nx,
       numeric::array& numeric_y, double ymin, double ymax, long ny,
       long num_entries, numeric::array& freq)
{
  double* x = reinterpret_cast<double *>
    (reinterpret_cast<PyArrayObject*>(numeric_x.ptr())->data);
  double* y = reinterpret_cast<double *>
    (reinterpret_cast<PyArrayObject*>(numeric_y.ptr())->data);
  double deltax = (xmax - xmin)/ nx / 2.0;
  double deltay = (ymax - ymin)/ ny / 2.0;

  long xbin, ybin;
  for(int i=0; i<num_entries; ++i) {
    xbin = static_cast<int>(ceil((x[i] - (xmin-deltax))/
				 (xmax-xmin)*nx*nx/(nx+1.0)));
    ybin = static_cast<int>(ceil((y[i] - (ymin-deltay))/
				 (ymax-ymin)*ny*ny/(ny+1.0)));
    if ((xbin>0) && (xbin<=nx) && (ybin>0) && (ybin<=ny)) {
      freq[make_tuple(ybin-1,xbin-1)] += 1;
    }
  }
}

BOOST_PYTHON_MODULE(hist2d)
{
  numeric::array::set_module_and_type("Numeric", "ArrayType");
  def("hist2d",&hist2d);
}

