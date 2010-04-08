#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/independent_operation.h"
#include "chef_elements_fixture.h"
#include "bunch_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Chef_elements_fixture)
{
    Chef_propagate_operation chef_propagate_operation(chef_elements);
}

BOOST_FIXTURE_TEST_CASE(apply, Chef_elements_fixture)
{
    Bunch_fixture b;

    Chef_propagate_operation chef_propagate_operation(chef_elements);
    //    multi_array_print(b.bunch.get_local_particles(), "particles before");
    chef_propagate_operation.apply(b.bunch);
    //    multi_array_print(b.bunch.get_local_particles(), "particles after");
}
// test_note: We need to check that apply actual produces the correct results.
//            As of this writing, it almost certainly doesn't
