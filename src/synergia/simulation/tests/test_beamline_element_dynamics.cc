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

#include "beamline/beamline.h"
#include "beamline/RefRegVisitor.h"

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(drift_propagation)
{
    const double pc = 2.784435311;
    const double drift_length = 1.0;
    const double x_offset = 0.1;
    const double dpp_offset = 0.04;

    // create chef lattice
    Lattice_sptr lattice_sptr(new Lattice("my_lattice"));

    Lattice_element my_drift("drift", "my_drift");
    my_drift.set_double_attribute("l", drift_length);
    lattice_sptr->append(my_drift);
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

    // calibrate the times
    RefRegVisitor registrar(proton);
    bml_ptr->accept(registrar);

    // test straight drift central particle
    bml_ptr->propagate(proton);
#if 0
    std::cout << "Proton state: " << proton.get_x() << " " << proton.get_y() << " " << proton.get_cdt() << " " << proton.get_npx() << " " << proton.get_npy() << " " << proton.get_ndp() << std::endl;
#endif

    BOOST_CHECK(floating_point_equal(proton.get_x(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_y(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_cdt(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npx(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npy(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_ndp(), 0.0, tolerance));

    // test drift from -xoffset to +xoffset
    proton.setStateToZero();

    const double proton_path_length = 
        std::sqrt(drift_length*drift_length + 4.0*x_offset*x_offset);
    const double xprime = 2.0*x_offset/proton_path_length;

    proton.set_x(-x_offset);
    proton.set_npx(xprime);
    bml_ptr->propagate(proton);

#if 0
    std::cout << "Proton state: " << proton.get_x() << " " << proton.get_y() << " " << proton.get_cdt() << " " << proton.get_npx() << " " << proton.get_npy() << " " << proton.get_ndp() << std::endl;
#endif

    double beta = reference_particle.get_beta();

    BOOST_CHECK(floating_point_equal(proton.get_x(), x_offset, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_y(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_cdt(),
                                     (proton_path_length/beta - drift_length/beta),
                                     tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npx(), xprime, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npy(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_ndp(), 0.0, tolerance));

    // check drift at offset dpp
    proton.setStateToZero();
    proton.set_ndp(dpp_offset);
    double proton_momentum = pc * (1.0 + dpp_offset);
    // make sure it agrees with CHEF's calculation
    BOOST_CHECK(floating_point_equal(proton_momentum, proton.Momentum(),
                                     tolerance));
    double proton_beta = proton.Beta();
    bml_ptr->propagate(proton);
    
#if 0
    std::cout << "Proton state: " << proton.get_x() << " " << proton.get_y() << " " << proton.get_cdt() << " " << proton.get_npx() << " " << proton.get_npy() << " " << proton.get_ndp() << std::endl;
#endif

    BOOST_CHECK(floating_point_equal(proton.get_x(),0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_y(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_cdt(),
                                     drift_length*(1.0/proton_beta-1.0/beta),
                                     tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npx(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_npy(), 0.0, tolerance));
    BOOST_CHECK(floating_point_equal(proton.get_ndp(), dpp_offset, tolerance));
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
}
