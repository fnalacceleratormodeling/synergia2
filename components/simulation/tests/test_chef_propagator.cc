#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "bunch_fixture.h"
#include "chef_elements_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;
const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Chef_elements_fixture)
{
    Chef_propagator chef_propagator(chef_elements);
}

BOOST_AUTO_TEST_CASE(apply)
{
    Chef_elements_fixture c;
    Bunch_fixture b;

    Chef_propagator chef_propagator(c.chef_elements);

    //    multi_array_print(bunch.get_local_particles(),"particles before");
    chef_propagator.apply(b.bunch);
    //    multi_array_print(bunch.get_local_particles(),"particles after");
}
// test_note: the apply test just verifies that the method doesn't crash.
//            At least one quantitative test would be desirable.
