#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/independent_operation.h"
#include "mapping_fixture.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Mapping_fixture)
{
    Fast_mapping fast_mapping(mapping, mapping_length);
    Fast_mapping_operation fast_mapping_operation(fast_mapping);
}

BOOST_FIXTURE_TEST_CASE(get_type, Mapping_fixture)
{
    Fast_mapping fast_mapping( mapping, mapping_length);
    Fast_mapping_operation fast_mapping_operation(fast_mapping);
    BOOST_CHECK_EQUAL(fast_mapping_operation.get_type(), "fast_mapping");
}

BOOST_FIXTURE_TEST_CASE(apply, Mapping_fixture)
{
    Fast_mapping fast_mapping( mapping, mapping_length);
    Fast_mapping_operation fast_mapping_operation(fast_mapping);
    //    multi_array_print(b.bunch.get_local_particles(), "particles before");
    const int verbosity = 5;
    Logger logger(0);
    fast_mapping_operation.apply(b.bunch, verbosity, logger);
    //    multi_array_print(b.bunch.get_local_particles(), "particles after");
}
// test_note: We need to check that apply actual produces the correct results.
//            As of this writing, it almost certainly doesn't

BOOST_FIXTURE_TEST_CASE(serialize_xml, Mapping_fixture)
{
    Fast_mapping fast_mapping(mapping, mapping_length);
    Fast_mapping_operation fast_mapping_operation(fast_mapping);

    xml_save<Fast_mapping_operation > (fast_mapping_operation,
            "fast_mapping_operation.xml");

    Fast_mapping_operation loaded;
    xml_load<Fast_mapping_operation > (loaded, "fast_mapping_operation.xml");
}

