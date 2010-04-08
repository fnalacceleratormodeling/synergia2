#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/independent_operation.h"
#include "mapping_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping);
    Fast_mapping_operation fast_mapping_operation(fast_mapping);
}

BOOST_FIXTURE_TEST_CASE(apply, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping);
    Fast_mapping_operation fast_mapping_operation(fast_mapping);
    //    multi_array_print(b.bunch.get_local_particles(), "particles before");
    fast_mapping.apply(b.bunch);
    //    multi_array_print(b.bunch.get_local_particles(), "particles after");
}
