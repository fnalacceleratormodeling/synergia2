#include "wrap_containers.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"


using namespace boost::python;

BOOST_PYTHON_MODULE(wrap_containers_python)
{

 def("take_vector_int", &take_vector_int);
 def("take_list_int", &take_list_int);
 def("take_vector_dd", &take_vector_dd);
 
}

