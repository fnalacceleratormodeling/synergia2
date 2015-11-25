#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "synergia/lattice/madx_adaptor_map.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

BOOST_AUTO_TEST_CASE(construct)
{
    MadX_adaptor_map madx_adaptor_map;
}

BOOST_AUTO_TEST_CASE(serialize_xml)
{
    MadX_adaptor_map madx_adaptor_map;
    xml_save<MadX_adaptor_map > (madx_adaptor_map,
            "madx_adaptor_map.xml");

    MadX_adaptor_map loaded;
    xml_load<MadX_adaptor_map > (loaded, "madx_adaptor_map.xml");
}
