#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "field_domain.h"
#include "array_nd/vector_helper.h"

BOOST_AUTO_TEST_CASE( construct_circular_field_domain )
{
    double length = 3.0;
    double radius = 4.0;
    std::vector<int> grid_shape(vector3(4,4,4));
    bool z_periodic = false;
    
    Cylindrical_field_domain cfd(length,radius,grid_shape,z_periodic);
}
