#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics_writer.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct)
{
    Multi_diagnostics_writer multi_writer;
}

BOOST_AUTO_TEST_CASE(append)
{
    Multi_diagnostics_writer multi_writer;
    Diagnostics_writer_sptr dummy_diagnostics_sptr(new Diagnostics_writer);
    multi_writer.append(dummy_diagnostics_sptr);
}

BOOST_AUTO_TEST_CASE(push_back)
{
    Multi_diagnostics_writer multi_writer;
    Diagnostics_writer_sptr dummy_diagnostics_sptr(new Diagnostics_writer);
    multi_writer.push_back(dummy_diagnostics_sptr);
}

BOOST_AUTO_TEST_CASE(iterate)
{
    Multi_diagnostics_writer multi_writer;
    Diagnostics_writer_sptr
            dummy_diagnostics_sptr1(new Diagnostics_writer);
    Diagnostics_writer_sptr
            dummy_diagnostics_sptr2(new Diagnostics_writer);
    Diagnostics_writer_sptr
            dummy_diagnostics_sptr3(new Diagnostics_writer);
    multi_writer.push_back(dummy_diagnostics_sptr1);
    multi_writer.push_back(dummy_diagnostics_sptr2);
    multi_writer.push_back(dummy_diagnostics_sptr3);
    int count(0);
    for (Multi_diagnostics_writer::iterator it = multi_writer.begin(); it
            != multi_writer.end(); ++it) {
        ++count;
    }
    BOOST_CHECK_EQUAL(count, 3);
}
