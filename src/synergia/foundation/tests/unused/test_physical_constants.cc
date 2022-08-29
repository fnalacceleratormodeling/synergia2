#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <stdexcept>

#include "synergia/foundation/physical_constants.h"
#include "basic_toolkit/PhysicsConstants.h"

using namespace pconstants;

const double tolerance = 4.0e-14;

BOOST_AUTO_TEST_CASE(physical_constants)
{
    // The values of the constants *really* should agree

    // check that the two names for particle masses are the same
    BOOST_CHECK_EQUAL(mp, proton_mass);
    BOOST_CHECK_EQUAL(me, electron_mass);
    BOOST_CHECK_EQUAL(mmu, muon_mass);

    // check that the synergia masses here agree with the CHEF masses
    BOOST_CHECK_EQUAL(proton_mass, PH_NORM_mp);
    BOOST_CHECK_EQUAL(electron_mass, PH_NORM_me);
    BOOST_CHECK_EQUAL(muon_mass, PH_NORM_mmu);

    // check other constants agree between synergia and CHEF
    BOOST_CHECK_EQUAL(c, PH_MKS_c);
    BOOST_CHECK_EQUAL(e, PH_MKS_e);

    // these are computed but they should be *really* close
    BOOST_CHECK_CLOSE(rp, PH_MKS_rp, tolerance);
    BOOST_CHECK_CLOSE(re, PH_MKS_re, tolerance);
}
