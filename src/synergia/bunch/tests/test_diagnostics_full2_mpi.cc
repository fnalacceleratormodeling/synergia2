#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

#include "bunch_sptr_fixture.h"

const double tolerance = 1.0e-11;

BOOST_FIXTURE_TEST_CASE(construct, Bunch_sptr_fixture)
{
    Diagnostics_full2 diagnostics("diagnostics_full2_mpi.h5");
}


BOOST_FIXTURE_TEST_CASE(write_, Bunch_sptr_fixture)
{
    Diagnostics_full2 diagnostics("diagnostics_full2_mpi.h5");
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();
}

BOOST_FIXTURE_TEST_CASE(write_several, Bunch_sptr_fixture)
{
    Diagnostics_full2 diagnostics("diagnostics_full2_mpi.h5");
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update_and_write();
    diagnostics.update_and_write();
    diagnostics.update_and_write();
    diagnostics.update_and_write();
}
