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
    Diagnostics_particles diagnostics("dummy.h5");
}


BOOST_FIXTURE_TEST_CASE(write_, Bunch_sptr_fixture)
{
    Diagnostics_particles diagnostics("diagnostics_particles_mpi1.h5");
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();
}

BOOST_FIXTURE_TEST_CASE(write_several, Bunch_sptr_fixture)
{
    Diagnostics_particles diagnostics("diagnostics_particles_mpi2.h5");
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();

    diagnostics.update();
    diagnostics.write();

    diagnostics.update();
    diagnostics.write();

    diagnostics.update();
    diagnostics.write();
}

BOOST_FIXTURE_TEST_CASE(write_min_max, Bunch_sptr_fixture)
{
    Diagnostics_particles diagnostics(
            "diagnostics_particles_mpi_35.h5", 3, 5);
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();
}
