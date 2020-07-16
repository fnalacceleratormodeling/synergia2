#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/reference_particle.h"

#include "synergia/bunch/bunch.h"

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/madx_adaptor_map.h"
#include "synergia/lattice/lattice_element_slice.h"

#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/bunch_simulator.h"

#include "synergia/utils/floating_point.h"

BOOST_AUTO_TEST_CASE(orbit_bump_timing)
{
    const double drlen = 1.0;
    const double vkval = 0.001;
                         
    const double tolerance = 1.0e-13;

    MadX_adaptor_map_sptr mxam(new MadX_adaptor_map);
    Lattice_sptr lattice_sptr(new Lattice("model", mxam));

    Lattice_element dr("drift", "dr");
    dr.set_double_attribute("l", drlen);

    Lattice_element vkp("vkicker", "vkp");
    vkp.set_double_attribute("kick", vkval);
    Lattice_element vkm("vkicker", "vkm");
    vkm.set_double_attribute("kick", -vkval);

    // create line with an orbit bump
    // drift straight
    lattice_sptr->append(dr);
    // kick up
    lattice_sptr->append(vkp);
    // drift diagonal up
    lattice_sptr->append(dr);
    // kick down to straighten out
    lattice_sptr->append(vkm);
    // drift straight at offset y
    lattice_sptr->append(dr);
    // kick down to return to origin
    lattice_sptr->append(vkm);
    // drift diagonally down
    lattice_sptr->append(dr);
    // kick up to straighten out
    lattice_sptr->append(vkp);
    // drift forward at origin
    lattice_sptr->append(dr);

    const double ke = 0.8;
    Reference_particle refpart(1, pconstants::mp, pconstants::mp + ke);
    double beta = refpart.get_beta();

    lattice_sptr->set_reference_particle(refpart);
    lattice_sptr->set_all_string_attribute("extractor_type", "libff");
    Commxx_sptr commxx(new Commxx);

    Bunch_sptr bunch_sptr(new Bunch(refpart, 1, 1.0e10, commxx));

    Bunch_simulator bunch_simulator(bunch_sptr);

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));

    // this sets the reference time of all lattice element slices
    MArray1d oneturn(stepper_sptr->get_lattice_simulator().tune_linear_lattice());

    Propagator propagator(stepper_sptr);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    propagator.propagate(bunch_simulator, 1, 1, 10);

    for (int i=0; i<4; ++i) {
        std::cout << std::setprecision(16) << "orbit propagation " << i << ", " << bunch_sptr->get_local_particles()[0][i] << std::endl;
        BOOST_CHECK(floating_point_equal(bunch_sptr->get_local_particles()[0][i], 0.0, co_tolerance));
    }

    const double py = vkval;
    const double pz = std::sqrt(1.0 - py*py);
    const double dy = (py/pz) * drlen;
    const double diag_length = std::sqrt(drlen*drlen + dy*dy);

    // calculate total pathlength
    // first drift straight
    double pathlen = drlen;
    // drift up diagonally
    pathlen += diag_length;
    // drift straight at offset y
    pathlen += drlen;
    // drift diagonally down
    pathlen += diag_length;
    // drift straight at 0 offset
    pathlen += drlen;

    double total_cdt = pathlen/beta;

    std::cout << std::setprecision(16) << "cdt from pathlength: " << total_cdt << std::endl;
    std::cout << std::setprecision(16) << "cdt from sum get_reftime: " << oneturn[Bunch::cdt] << std::endl;

    BOOST_CHECK(floating_point_equal(total_cdt, oneturn[Bunch::cdt], 1.0e-12));
}
