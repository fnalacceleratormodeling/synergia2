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

#include "beamline/beamline.h"
#include "beamline/RefRegVisitor.h"

#include "synergia/utils/boost_test_mpi_fixture.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture) // needed to initialize MPI

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(drift_propagation)
{
    const double pc = 2.784435311;
    const double drift_length = 1.0;
    const double real_particles = 1.0e9;
    const int map_order = 1;    // not really used for step for test

    const double x_offset = 0.1;
    const double dpp_offset = 0.04;

    // create chef lattice
    Lattice_sptr lattice_sptr(new Lattice("my_lattice"));

    Lattice_element my_drift("drift", "my_drift");
    my_drift.set_double_attribute("l", drift_length);
    my_drift.set_string_attribute("extractor_type", "libff");
    // my_drift.set_string_attribute("extractor_type", "chef_propagate");
    lattice_sptr->append(my_drift);
#if 1
    std::cout << "this is my lattice" << std::endl;
    std::cout << lattice_sptr->as_string() << std::endl;
    std::cout << "end of my lattice" << std::endl;
#endif

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);
    std::cout << "before commxx" << std::endl;
    Commxx_sptr commxx(new Commxx());
    std::cout << "after ocmmxx" << std::endl;
    // 3 particles in the test + reference
    Bunch_sptr bunch_sptr(new Bunch(reference_particle, 4, real_particles, commxx));
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
    local_particles[1][Bunch::x] = 0.0;
    local_particles[1][Bunch::xp] = 0.0;
    local_particles[1][Bunch::y] = 0.0;
    local_particles[1][Bunch::yp] = 0.0;
    local_particles[1][Bunch::cdt] = 0.0;
    local_particles[1][Bunch::dpop] = 0.0;

    // test drift from -xoffset to +xoffset
    const double proton_path_length =
        std::sqrt(drift_length*drift_length + 4.0*x_offset*x_offset);
    const double xprime = 2.0*x_offset/proton_path_length;
    local_particles[2][Bunch::x] = (-x_offset);
    local_particles[2][Bunch::xp] = (xprime);
    local_particles[2][Bunch::y] = 0.0;
    local_particles[2][Bunch::yp] = 0.0;
    local_particles[2][Bunch::cdt] = 0.0;
    local_particles[2][Bunch::dpop] = 0.0;

    // check drift at offset dpp
    local_particles[3][Bunch::x] = 0.0;
    local_particles[3][Bunch::xp] = 0.0;
    local_particles[3][Bunch::y] = 0.0;
    local_particles[3][Bunch::yp] = 0.0;
    local_particles[3][Bunch::cdt] = 0.0;
    local_particles[3][Bunch::dpop] = dpp_offset;

    propagator.propagate(bunch_simulator, 1);

    // subtract reference particle cdt
    for (int j=1; j<4; ++j) {
        local_particles[j][Bunch::cdt] -= local_particles[0][Bunch::cdt];
    }
    std::cout << "propagating 0 particle through drift" << std::endl;
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::x], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::xp], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::y], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::yp], 0.0, tolerance));
    std::cout << "particle cdt (should be 0): " << local_particles[1][Bunch::cdt] << std::endl;
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::cdt], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[1][Bunch::dpop], 0.0, tolerance));

    // no momentum offset so particle beta is the same as the reference particle beta
    double beta = reference_particle.get_beta();

    std::cout << "propagating particle at slat through drift" << std::endl;
    std::cout << "particle x (should be " << x_offset  << "): " << local_particles[2][Bunch::x] << std::endl;
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::x], x_offset, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::y], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::cdt],
                                     (proton_path_length/beta - drift_length/beta),
                                     tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::xp], xprime, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::yp], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[2][Bunch::dpop], 0.0, tolerance));

    // particle has momentum offset so calculate its beta from its momentum
    double proton_momentum = pc * (1.0 + dpp_offset);
    Four_momentum test_four_momentum(pconstants::proton_mass);
    test_four_momentum.set_momentum(proton_momentum);
    double proton_beta = test_four_momentum.get_beta();

    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::x],0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::y], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::cdt],
                                     drift_length*(1.0/proton_beta-1.0/beta),
                                     tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::xp], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::yp], 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(local_particles[3][Bunch::dpop], dpp_offset, tolerance));
}

BOOST_AUTO_TEST_CASE(sbend_propagation)
{
    const double pc = 2.784435311;
    const double sbend_length = 1.09812;
    const double sbend_angle = 0.01567741327;
    const double dpp_offset = 0.04;
    const double sbend_tolerance = 1.0e-7;

    // create chef lattice
    Lattice_sptr lattice_sptr(new Lattice("my_lattice"));

    Lattice_element my_sbend("sbend", "my_sbend");
    my_sbend.set_double_attribute("l", sbend_length);
    my_sbend.set_double_attribute("angle", sbend_angle);
    lattice_sptr->append(my_sbend);
#if 0
    lattice_sptr->print();
#endif

    Four_momentum four_momentum(pconstants::proton_mass);
    four_momentum.set_momentum(pc);

    Reference_particle reference_particle(1, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

    Lattice_simulator lattice_simulator(lattice_sptr, 1);
    BmlPtr bml_ptr(lattice_simulator.get_chef_lattice().get_beamline_sptr());

    double total_energy = reference_particle.get_total_energy();
    Proton proton(total_energy);
    double e0 = total_energy;
    double p0 = pc;
    double beta = proton.Beta();

    // calibrate the times
    RefRegVisitor registrar(proton);
    bml_ptr->accept(registrar);

    // the central particle should just go through the sbend
    bml_ptr->propagate(proton);
    BOOST_CHECK(floating_point_equal(proton.get_x(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_y(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_cdt(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npx(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npy(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_ndp(), 0.0, tolerance));

    // calculate where it should be

    double brho = p0 * 1.0e9/pconstants::c;
    double magnetic_field = sbend_angle*brho/sbend_length;

    double radius = brho/magnetic_field; // radius of curvature for on-momentum

    /*
     * on-momentum particle moves around the circle of radius
     * Make the center of the circle (0, -radius)
     * the nominal particle travels length l.  The sbend edge will be
     * the line that includes the (0, -radius) and the end of travel.
     */
    double sbend_exit_x = radius * std::sin(sbend_angle);
    double sbend_exit_y = -radius * (1.0 - std::cos(sbend_angle));
    
    // check propagation of off momentum particle
    proton.setStateToZero();
    proton.set_ndp(dpp_offset);
    double beta1 = proton.Beta();
    double p1 = p0 * (1.0+dpp_offset);
    // check momentum that it agrees with CHEF momentum
    BOOST_CHECK(floating_point_equal(p1, proton.Momentum(),
                                     tolerance));

    bml_ptr->propagate(proton);
    double radius1 = (p1/magnetic_field) * (1.0e9/pconstants::c);
    double delta_radius = radius1-radius;
    double cot_theta = std::cos(sbend_angle)/std::sin(sbend_angle);
    /*
     * particle 1 exit point at the intersection of the circle of radius
     * radius1 centered at
     * (0, -radius1) with line from (0,-radius), (sbend_exit_x, sbend_exit_y)
     * which works out to be
     * y = -radius + cos(angle)/sin(angle)*x
     *
     *  intersection with circle
     * 
     * x**2 + (y + radius1)**2 = radius1**2
     *
     * substitute
     * (1+ct**2)*x**2 + 2*ct*dr*x + dr**2-radius1**2 = 0 and
     * solve quadratic equation
     *
     * px/pref calculation L. Michelotti:
     * px/pref = (dp/pref) * sin(theta_0)
     *
     */
    double a = 1.0 + cot_theta*cot_theta;
    double b = 2.0*cot_theta*delta_radius;
    double c = delta_radius*delta_radius - radius1*radius1;
    // quadratic solution (for 1/x, inverted)
    double sbend_exit_x1 = 2*c/(-b - std::sqrt(b*b - 4*a*c));
    double sbend_exit_y1 = -radius + cot_theta*sbend_exit_x1;
    double d_offset = std::sqrt(std::pow(sbend_exit_x-sbend_exit_x1,2) +
                                std::pow(sbend_exit_y-sbend_exit_y1,2));
    // offset distance is the x offset in accelerator coordinates
    BOOST_CHECK(floating_point_equal(d_offset, proton.get_x(), sbend_tolerance));
    // check px/pref
    double px_over_pref = dpp_offset * std::sin(sbend_angle);
    BOOST_CHECK(floating_point_equal(proton.get_npx(), px_over_pref, sbend_tolerance));
    // check path length difference
    // path length for particle 0 is just sbend_length
    // calculate path length for particle 1 on its circular trajectory
    // origin of circle is (0, -radius1)
    // vector of particle start position is (0, radius1)
    // vector of particle end position is (sbend_exit_x1, sbend_exit_y1-(-radius1))
    // calculate dot product to get angle
    double v0dotv1 = sbend_exit_y1+radius1; // *radius1 which will cancel
    double v1_norm = std::sqrt(pow(sbend_exit_x1, 2) +
                               pow(sbend_exit_y1+radius1, 2));
    double theta1 = std::acos(v0dotv1/v1_norm);  // * 1/radius1
    double sbend_length1 = theta1 * radius1;
    double delta_cdt = sbend_length1/beta1 - sbend_length/beta;
    BOOST_CHECK(floating_point_equal(proton.get_cdt(), delta_cdt, sbend_tolerance));

    // check propagation of x-offset particle on-momentum.
    proton.setStateToZero();
    double x_offset = .001;
    proton.set_x(x_offset);
    bml_ptr->propagate(proton);

    // In this case, the trajectory follows a circle of original radius
    // but centered on (0, -R+x_offset).  The exit point is the intersection of that
    // circle with the line of close cot(theta) passing through (0, -R) which is the solution
    // of the quadratic equation
    // (1+ct**2)*x**2 - 2*ct*x_offset*x + x_offset**2-radius**2 = 0
    a = 1.0 + cot_theta*cot_theta;
    b = -2.0*cot_theta*x_offset;
    c = x_offset*x_offset - radius*radius;
    // quadratic solution (for 1/x, inverted)
    double sbend_exit_x2 = 2*c/(-b - std::sqrt(b*b - 4*a*c));
    double sbend_exit_y2 = -radius + cot_theta*sbend_exit_x2;
    d_offset = std::sqrt(std::pow(sbend_exit_x-sbend_exit_x2,2) +
                                std::pow(sbend_exit_y-sbend_exit_y2,2));
    // offset distance is the x offset in accelerator coordinates
#if 0
    std::cout << "central  exit position: " << sbend_exit_x  << ", " << sbend_exit_y  << std::endl;
    std::cout << "particle exit position: " << sbend_exit_x2 << ", " << sbend_exit_y2 << std::endl;
    std::cout << "proton.get_x(): " << proton.get_x() << ", d_offset: " << d_offset << std::endl;
#endif
    BOOST_CHECK(floating_point_equal(proton.get_x(), d_offset, sbend_tolerance));

    // get pathlength
    // start position vector is (0, -radius)
    // final position vector is (sbend_exit_x2, sbend_exit_y2-(-radius+x_offset)
    double v0dotv2 = sbend_exit_y2+radius-x_offset; // times radius
    double v2norm = std::sqrt(pow(sbend_exit_x2, 2) +
    		                  pow(sbend_exit_y2-radius+x_offset, 2));
    double theta2 = std::acos(v0dotv2/radius);
    double sbend_length2 = theta2 * radius;
    delta_cdt = sbend_length2/beta - sbend_length/beta;
#if 0
    std::cout << "proton.get_cdt(): " << proton.get_cdt() << ", delta_cdt: " << delta_cdt << std::endl;
#endif
    BOOST_CHECK(floating_point_equal(proton.get_cdt(), delta_cdt, sbend_tolerance));

    // check x'.  Transverse momentum is the sin of the angle between the
    // vector normal to the sbend exit face and the particle trajectory vector.
    // calculate with the cross product.
    double exit_normal_x = (sbend_exit_y2-(-radius+x_offset))/radius;
    double exit_normal_y = -sbend_exit_x2/radius;
    double face_normal_x = std::cos(sbend_angle);
    double face_normal_y = -std::sin(sbend_angle);
    // cross product
    double exit_px = face_normal_x*exit_normal_y - face_normal_y*exit_normal_x;
#if 0
    std::cout << "exit_px: " << exit_px << ", proton.get_npx(): " << proton.get_npx() << std::endl;
#endif
    BOOST_CHECK(floating_point_equal(proton.get_npx(), exit_px, sbend_tolerance));
}
