#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/floating_point.h"
#include "synergia/utils/serialization.h"

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/propagator.h"

#include "synergia/bunch/bunch.h"

#include "synergia/utils/boost_test_mpi_fixture.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture); // needed to initialize MPI

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(drift_propagation)
{
    const double pc = 2.784435311;
    const double drift_length = 1.0;
    const double real_particles = 1.0e9;
    const int map_order = 1;    // not really used for step for test

    const double x_offset = 0.1;
    const double dpp_offset = 0.04;

    const double qoff = 0.01;
    const double vrel = 0.01;

    // create chef lattice
    Lattice_sptr lattice_sptr(new Lattice("my_lattice"));

    Lattice_element my_drift("drift", "my_drift");
    my_drift.set_double_attribute("l", drift_length);
    my_drift.set_string_attribute("extractor_type", "libff");
    // my_drift.set_string_attribute("extractor_type", "chef_propagate");
    lattice_sptr->append(my_drift);

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    Commxx_sptr commxx(new Commxx());

    // 3 particles in the test + reference
    // pad to 8 particles to satisfy possible AVX512 alignment requirement
    Bunch_sptr bunch_sptr(new Bunch(reference_particle, 8, real_particles, commxx));
    Bunch_simulator bunch_simulator(bunch_sptr);

    MArray2d_ref local_particles(bunch_sptr->get_local_particles());

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, map_order, 1));
    Propagator propagator(stepper_sptr);

    // reference particle at 0
    MArray1d ref_state(reference_particle.get_state());
    for (int i=0; i<6; ++i) {
        local_particles[0][i] = ref_state[i];
    }

    // test straight drift central particle
    local_particles[1][Bunch::x]    = 0.0;
    local_particles[1][Bunch::xp]   = 0.0;
    local_particles[1][Bunch::y]    = 0.0;
    local_particles[1][Bunch::yp]   = 0.0;
    local_particles[1][Bunch::cdt]  = 0.0;
    local_particles[1][Bunch::dpop] = 0.0;

    // initial values
    const double pr = bunch_sptr->get_reference_particle().get_momentum();
    const double m  = bunch_sptr->get_mass();
    const double q  = qoff;
    const double l  = drift_length;
    const double g  = 1.0 / ( 1.0 - vrel );
    const double g2 = 1.0 / ( (1.0-vrel) * (1.0-vrel) );
    const double beta = reference_particle.get_beta();

    const double xp2 = q * m / sqrt( l*l*m*m - q*q*pr*pr );
    const double xp3 = q * m / sqrt( l*l*m*m - 2*q*q*pr*pr );
    const double xp4 = q * m * g * sqrt( (l*l*(1-g2)*pr*pr + m*m) / (l*l*(1-g2)*pr*pr+m*m) / (l*l*(1-g2)*pr*pr+m*m-q*q*g*g*pr*pr) );

    const double cdt4 = vrel * l / (2*beta);

    const double dp2 = sqrt( xp2*xp2*(pr*pr+m*m)/(m*m) + 1 ) - 1;
    const double dp3 = sqrt( 2*xp3*xp3*(pr*pr+m*m)/(m*m) + 1 ) - 1;
    const double dp4 = sqrt( (xp4*xp4*(pr*pr+m*m) + g*g*m*m)/((1-g*g)*pr*pr+m*m) ) - 1;

    //std::cout << "xp2 = " << xp2 << ", xp3 = " << xp3 << ", xp4 = " << xp4 << "\n";
    //std::cout << "cdt4 = " << cdt4 << "\n";
    //std::cout << "dp2 = " << dp2 << ", dp3 = " << dp3 << ", dp4 = " << dp4 << "\n";

    local_particles[2][Bunch::x]    = -q;
    local_particles[2][Bunch::xp]   = xp2;
    local_particles[2][Bunch::y]    = 0.0;
    local_particles[2][Bunch::yp]   = 0.0;
    local_particles[2][Bunch::cdt]  = 0.0;
    local_particles[2][Bunch::dpop] = dp2;

    // check drift at offset dpp
    local_particles[3][Bunch::x]    = -q;
    local_particles[3][Bunch::xp]   = xp3;
    local_particles[3][Bunch::y]    = 0.0;
    local_particles[3][Bunch::yp]   = xp3;
    local_particles[3][Bunch::cdt]  = 0.0;
    local_particles[3][Bunch::dpop] = dp3;

    local_particles[4][Bunch::x]    = -q/2.0;
    local_particles[4][Bunch::xp]   = xp4;
    local_particles[4][Bunch::y]    = 0.0;
    local_particles[4][Bunch::yp]   = 0.0;
    local_particles[4][Bunch::cdt]  = cdt4;
    local_particles[4][Bunch::dpop] = dp4;

#if 0
    std::cout << std::setprecision(10);

    for (int i=0; i<5; ++i)
    {
        std::cout << "particle " << i << ": ";
        for(int j=0; j<6; ++j)
        {
            std::cout << local_particles[i][j] << ", ";
        }
        std::cout << "\n";
    }
#endif


    propagator.propagate(bunch_simulator, 1);

#if 0
    std::cout << std::setprecision(10);

    for (int i=0; i<5; ++i)
    {
        std::cout << "particle " << i << ": ";
        for(int j=0; j<6; ++j)
        {
            std::cout << local_particles[i][j] << ", ";
        }
        std::cout << "\n";
    }
#endif

    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::x],    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::xp],   0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::y],    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::yp],   0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::cdt],  0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::dpop], 0.0, tolerance));

    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::x],    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::xp],   xp2, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::y],    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::yp],   0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::cdt],  0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::dpop], dp2, tolerance));

    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::x],    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::xp],   xp3, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::y],    qoff, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::yp],   xp3, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::cdt],  0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::dpop], dp3, tolerance));

    BOOST_CHECK(floating_point_equal(local_particles[4][Bunch::x],    q/2.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[4][Bunch::xp],   xp4,   tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[4][Bunch::y],    0.0,   tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[4][Bunch::yp],   0.0,   tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[4][Bunch::cdt],  -cdt4, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[4][Bunch::dpop], dp4,   tolerance));
}


