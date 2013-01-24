#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "synergia/lattice/mad8_adaptor_map.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct)
{
    Mad8_adaptor_map mad8_adaptor_map;
}

BOOST_AUTO_TEST_CASE(serialize_xml)
{
    Mad8_adaptor_map mad8_adaptor_map;
    xml_save<Mad8_adaptor_map > (mad8_adaptor_map,
            "mad8_adaptor_map.xml");

    Mad8_adaptor_map loaded;
    xml_load<Mad8_adaptor_map > (loaded, "mad8_adaptor_map.xml");
}
