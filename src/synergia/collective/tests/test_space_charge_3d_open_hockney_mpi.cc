#include "synergia/utils/catch.hpp"

#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/foundation/physical_constants.h"

TEST_CASE("space charge", "[Space_charge_3d_open_hockney]")
{
    CHECK(1);
}


const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.7e11;
const int total_num = 10000;
const double total_energy = 125.0;

// 1 + (odd number) * 8 + 2 to expand the domein on both sides
const int rod_num_particles = 1 + 250001*8; 

const double rod_lowgamma = 61.0/60.0;
const double rod_highgamma = 61.0/11.0;
const double rod_real_num = 5.0e9;
const double rod_length = 0.1;

// radius of rod particles
const double rod_radius = 1.0e-6; 

// radius of the probe particle
const double rod_probe = 1.0e-3;  

struct Rod_bunch_fixture_lowgamma
{
    Rod_bunch_fixture_lowgamma() 
        : bsim(Bunch_simulator::create_single_bunch_simulator(
                Reference_particle(charge, Four_momentum(mass, mass*rod_lowgamma)),
                rod_num_particles, rod_real_num))
    {
        auto& bunch = bsim.get_bunch();

        //bunch.set_longitudinal_boundary(LongitudinalBoundary::periodic, rod_length);
        auto local_particles = bunch.get_local_particles();

        // a ring of 8 particles around each longitudinal location
        int num_longitudinal = (rod_num_particles-1)/8;
        double dz = rod_length/(num_longitudinal-1);
        double r2o2 = std::sqrt(2.0)/2.0;
        double z = -rod_length/2.0;

        auto const& ref = bunch.get_reference_particle();
        double rod_beta = ref.get_beta();

        for (int i=1; i<rod_num_particles; i+=8, z+=dz) 
        {
            local_particles(i, Bunch::x) = rod_radius;
            local_particles(i, Bunch::y) = 0.0;

            local_particles(i+1, Bunch::x) = rod_radius*r2o2;
            local_particles(i+1, Bunch::y) = rod_radius*r2o2;

            local_particles(i+2, Bunch::x) = 0.0;
            local_particles(i+2, Bunch::y) = rod_radius;

            local_particles(i+3, Bunch::x) = -rod_radius*r2o2;
            local_particles(i+3, Bunch::y) =  rod_radius*r2o2;

            local_particles(i+4, Bunch::x) = -rod_radius;
            local_particles(i+4, Bunch::y) = 0.0;

            local_particles(i+5, Bunch::x) = -rod_radius*r2o2;
            local_particles(i+5, Bunch::y) = -rod_radius*r2o2;

            local_particles(i+6, Bunch::x) = 0.0;
            local_particles(i+6, Bunch::y) = -rod_radius;

            local_particles(i+7, Bunch::x) = rod_radius*r2o2;
            local_particles(i+7, Bunch::y) = -rod_radius*r2o2;


            for (int j=i; j<i+8; ++j) 
            {
                local_particles(j, Bunch::cdt) = z/rod_beta;
                local_particles(j, Bunch::xp) = 0.0;
                local_particles(j, Bunch::yp) = 0.0;
                local_particles(j, Bunch::dpop) = 0.0;
                local_particles(j, Bunch::id) = j;
            }
        }

        // when the probe is too far, it falls outside the domain and does
        // not get any sc kicks.
        local_particles(0, Bunch::x) = 80*rod_radius;
        local_particles(0, Bunch::y) = 0.0;
        local_particles(0, Bunch::xp) = 0.0;
        local_particles(0, Bunch::yp) = 0.0;
        local_particles(0, Bunch::cdt) = 0.0;
        local_particles(0, Bunch::dpop) = 0.0;
        local_particles(0, Bunch::id) = 0.0;

        // check in particles
        bunch.checkin_particles();
    }

    Bunch_simulator bsim;
};



TEST_CASE("real_apply_full_lowgamma", "[Rod_bunch]")
{
    auto logger = Logger(0, LoggerV::DEBUG);
    auto simlogger = Logger(0, LoggerV::INFO_STEP);

    const int gridx = 256;
    const int gridy = 256;
    const int gridz = 64;

    const double time_fraction = 1.0;
    const double step_length = 0.1;

    // bunch
    Rod_bunch_fixture_lowgamma fixture;

    auto&       bunch = fixture.bsim.get_bunch();
    auto const& ref   = bunch.get_reference_particle();
    auto        parts = bunch.get_local_particles();

    const double beta = ref.get_beta();
    const double gamma = ref.get_gamma();
    const double betagamma = beta * gamma;

    const double time_step = step_length/(beta*pconstants::c);
    //const double bunchlen = bunch.get_longitudinal_boundary().second;
    const double bunchlen = 0.1;

    // check out particles before print
    bunch.checkout_particles();

    // print intital coordinates
    logger << "real_apply_full_lowgamma first four particles (x y z):" << std::endl;
    for (int k=0; k<4; ++k) 
    {
        logger << k << ": " 
            << parts(k, 0) << ", " 
            << parts(k, 2) << ", " 
            << parts(k, 4) << std::endl;
    }

    logger << "last four particles (x y z):" << std::endl;
    for (int k=bunch.get_local_num()-4; k<bunch.get_local_num(); ++k) 
    {
        logger << k << ": " 
            << parts(k, 0) << ", " 
            << parts(k, 2) << ", " 
            << parts(k, 4) << std::endl;
    }

    logger << std::endl;

    // space charge operator
    auto sc_ops = Space_charge_3d_open_hockney_options(gridx, gridy, gridz);
    sc_ops.comm_group_size = 1;
    sc_ops.green_fn = green_fn_t::linear;

    auto sc = Space_charge_3d_open_hockney(sc_ops);

    // set domain
    std::array<double, 3> offset = {0, 0, 0};
    std::array<double, 3> size = { parts(0, 0) * 4, parts(0, 0)*4, bunchlen/beta };
    sc.set_fixed_domain(offset, size);

    // apply space charge operator
    sc.apply(fixture.bsim, time_step, simlogger);

    // check out particles
    bunch.checkout_particles();

    // print
    logger << "after sc::apply : bunch.local_particles(0, 0): " << parts(0, 0) << std::endl;

    // Rod of charge Q over length L
    // E field at radius r $$ E = \frac{1}{2 \pi \epsilon_0} \frac{Q}{L} \frac{1}{r} $$
    // B field at radius r $$ E = \frac{\mu_0}{2 \pi } \frac{Q v}{L} \frac{1}{r} $$
    // Net EM force on electric+magnetic on probe of charge q from E-B cancellation
    // $$ F = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{1}{r}
    // travel over distance D at velocity v
    // \frac{\Delta p}{p} = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{D}{m v^2} \frac{1}{r}
    // convert to usual units
    // \frac{\Delta p}{p} = \frac{2 N r_p}{L \beta^2 \gamma^3} \frac{D}{r}

    double L = bunchlen;
    double N = bunch.get_real_num();

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "beta: " << beta << std::endl;
    logger << "gamma: " << gamma << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << parts(0, Bunch::x) << std::endl;

    double computed_dpop = ((2.0*N*pconstants::rp)/(L*betagamma*betagamma*gamma)) *
            (step_length/parts(0, Bunch::x));

    logger << "computed dpop: " << computed_dpop << std::endl;
    logger << "particle dpop: " << parts(0, 1) << std::endl;

    CHECK(parts(0, Bunch::xp) == Approx(computed_dpop).margin(.01));

    int nkicks = 0;
    for (int k=0; k<bunch.get_local_num(); ++k) 
    {
        if ((parts(k, 1)!=0.0) || (parts(k, 3)!=0.0)) 
        {
            ++nkicks;
            
            if (nkicks < 10) 
            {
                logger << "kick: " << nkicks 
                    << ", particle " << k << ": " 
                    << parts(k, 0) << ", " 
                    << parts(k, 1) << ", " 
                    << parts(k, 2) << ", " 
                    << parts(k, 3) << ", " 
                    << parts(k, 4) << ", " 
                    << parts(k, 5) << std::endl;
            }
        }
    }



#if 0
    Bunch original_bunch(bunch);

    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    bool show_output(false);
    Logger logger(0, show_output);

    logger << "real_apply_full_lowgamma first four particles (x y z):" << std::endl;
    for (int k=0; k<4; ++k) {
        logger << k<< ": " << bunch.get_local_particles()[k][0] << ", " <<
                bunch.get_local_particles()[k][2] << ", " <<
                bunch.get_local_particles()[k][4] << std::endl;
    }
    logger << "last four particles (x y z):" << std::endl;
    for (int k=bunch.get_local_num()-4; k<bunch.get_local_num(); ++k) {
        logger << k<<": " << bunch.get_local_particles()[k][0] << ", " <<
                bunch.get_local_particles()[k][2] << ", " <<
                bunch.get_local_particles()[k][4] << std::endl;
    }
    logger << std::endl;

    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, true, false, 0.0, false, 16.0);

    space_charge.update_domain(bunch);
    Rectangular_grid_domain_sptr orig_domain_sptr(space_charge.get_domain_sptr());
    std::vector<int> sc_grid_shape(orig_domain_sptr->get_grid_shape());
    std::vector<double> sc_size(orig_domain_sptr->get_physical_size());
    std::vector<double> sc_offs(orig_domain_sptr->get_physical_offset());

    std::vector<double> domain_sizezyx(3);
    std::vector<double> domain_offsetzyx(3);
    domain_offsetzyx[2] = 0.0;
    domain_offsetzyx[1] = 0.0;
    domain_offsetzyx[0] = 0.0;
    domain_sizezyx[2] = bunch.get_local_particles()[0][Bunch::x] * 4;
    domain_sizezyx[1] = domain_sizezyx[2];
    domain_sizezyx[0] = bunchlen/beta;
    std::vector<int> grid_shapezyx(3);
    grid_shapezyx[0] = grid_shape[2];
    grid_shapezyx[1] = grid_shape[1];
    grid_shapezyx[2] = grid_shape[0];

    Rectangular_grid_domain_sptr fixed_domain(
            new Rectangular_grid_domain(domain_sizezyx, domain_offsetzyx, grid_shapezyx));
    space_charge.set_fixed_domain(fixed_domain);

    Rectangular_grid_sptr local_charge_density(
            space_charge.get_local_charge_density(bunch));
    std::vector<int> local_rho_shape(local_charge_density->get_domain_sptr()->get_grid_shape());
    std::vector<double> local_rho_phys_size(local_charge_density->get_domain_sptr()->get_physical_size());
    std::vector<double> local_rho_phys_off(local_charge_density->get_domain_sptr()->get_physical_offset());
    std::vector<double> local_rho_left(local_charge_density->get_domain_sptr()->get_left());

    Step dummy_step(step_length);

    const int verbosity = 99;

    space_charge.apply(bunch, time_step, dummy_step, verbosity, logger);
    logger << "bunch final state: " << bunch.get_state() << std::endl;
    bunch.convert_to_state(Bunch::fixed_z_lab);
    logger << "after sc::apply : bunch.get_local_particles()[0][0]: " <<
            bunch.get_local_particles()[0][0] << std::endl;
    // Rod of charge Q over length L
    // E field at radius r $$ E = \frac{1}{2 \pi \epsilon_0} \frac{Q}{L} \frac{1}{r} $$
    // B field at radius r $$ E = \frac{\mu_0}{2 \pi } \frac{Q v}{L} \frac{1}{r} $$
    // Net EM force on electric+magnetic on probe of charge q from E-B cancellation
    // $$ F = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{1}{r}
    // travel over distance D at velocity v
    // \frac{\Delta p}{p} = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{D}{m v^2} \frac{1}{r}
    // convert to usual units
    // \frac{\Delta p}{p} = \frac{2 N r_p}{L \beta^2 \gamma^3} \frac{D}{r}

    double L = bunch.get_z_period_length();
    double N = bunch.get_real_num();
    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(L*betagamma*betagamma*gamma)) *
            (step_length/bunch.get_local_particles()[0][Bunch::x]);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);

    int nkicks = 0;
    for (int k=0; k<bunch.get_local_num(); ++k) {
        if ((bunch.get_local_particles()[k][1] != 0.0) ||
                (bunch.get_local_particles()[k][3] != 0.0)) {
            ++nkicks;
            if (nkicks < 10) {
                logger << "kick: " << nkicks << "particle " << k << ": " << bunch.get_local_particles()[k][0] << ", " <<
                          bunch.get_local_particles()[k][1] << ", " <<
                          bunch.get_local_particles()[k][2] << ", " <<
                          bunch.get_local_particles()[k][3] << ", " <<
                          bunch.get_local_particles()[k][4] << ", " <<
                          bunch.get_local_particles()[k][5] << std::endl;
            }
        }
    }
#endif

}



#if 0
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/floating_point.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/hdf5_file.h"
#include "gaussian_charge_density.h"
#include "space_charge_bunch_fixtures.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct1)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    Commxx_sptr comm_sptr(new Commxx);

    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    Commxx_sptr comm_sptr(new Commxx);
    bool longitudinal_kicks(false);
    bool periodic_z(true);
    double z_period(1.1);
    bool grid_entire_period(true);
    double n_sigma(7.0);

    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape,
            longitudinal_kicks, periodic_z, z_period, grid_entire_period,
            n_sigma);
}

BOOST_AUTO_TEST_CASE(construct_bad_period)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    Commxx_sptr comm_sptr(new Commxx);
    bool longitudinal_kicks(false);
    bool periodic_z(true);
    double z_period(0.0);
    bool grid_entire_period(true);
    double n_sigma(7.0);

    bool caught_error(false);
    try {
        Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape,
                longitudinal_kicks, periodic_z, z_period, grid_entire_period,
                n_sigma);
    }
    catch (std::runtime_error &) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);
}

BOOST_AUTO_TEST_CASE(get_n_sigma)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    Commxx_sptr comm_sptr(new Commxx);
    bool longitudinal_kicks(false);
    bool periodic_z(true);
    double z_period(1.1);
    bool grid_entire_period(true);
    double n_sigma(7.0);

    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape,
            longitudinal_kicks, periodic_z, z_period, grid_entire_period,
            n_sigma);
    BOOST_CHECK_CLOSE(space_charge.get_n_sigma(), n_sigma, tolerance);
}

BOOST_FIXTURE_TEST_CASE(update_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    space_charge.update_domain(bunch);
}

BOOST_FIXTURE_TEST_CASE(get_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    space_charge.update_domain(bunch);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[0],
            space_charge.get_n_sigma() * stdz, tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[1],
            space_charge.get_n_sigma() * stdy, tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[2],
            space_charge.get_n_sigma() * stdx, tolerance);
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[0],
                    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[1],
                    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[2],
                    0.0, tolerance));
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[0],
            grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[1],
            grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[2],
            grid_shape[0]);
}

BOOST_FIXTURE_TEST_CASE(get_doubled_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    space_charge.update_domain(bunch);
    BOOST_CHECK_EQUAL(space_charge.get_doubled_domain_sptr()->get_grid_shape()[0],
            2*grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_doubled_domain_sptr()->get_grid_shape()[1],
            2*grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_doubled_domain_sptr()->get_grid_shape()[2],
            2*grid_shape[0]);
}

BOOST_FIXTURE_TEST_CASE(set_fixed_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    space_charge.update_domain(bunch);
    Rectangular_grid_domain_sptr domain_sptr = space_charge.get_domain_sptr();
    double scale = 1.23;
    double shift = 0.2345;
    std::vector<double > scaled_size(3), shifted_offset(3);
    for (int i = 0; i < 3; ++i) {
        scaled_size[i] = scale * domain_sptr->get_physical_size()[i];
        shifted_offset[i] = shift + domain_sptr->get_physical_offset()[i];
    }
    Rectangular_grid_domain_sptr new_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    scaled_size, shifted_offset, domain_sptr->get_grid_shape(),
                    domain_sptr->is_periodic()));
    space_charge.set_fixed_domain(new_domain_sptr);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[0],
            scaled_size[0], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[1],
            scaled_size[1], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[2],
            scaled_size[2], tolerance);
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[0],
                    shifted_offset[0], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[1],
                    shifted_offset[1], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[2],
                    shifted_offset[2], tolerance));
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[0],
            grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[1],
            grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[2],
            grid_shape[0]);

    // make sure that update domain has no effect
    space_charge.update_domain(bunch);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[0],
            scaled_size[0], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[1],
            scaled_size[1], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain().get_physical_size()[2],
            scaled_size[2], tolerance);
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[0],
                    shifted_offset[0], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[1],
                    shifted_offset[1], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain().get_physical_offset()[2],
                    shifted_offset[2], tolerance));
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[0],
            grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[1],
            grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_domain().get_grid_shape()[2],
            grid_shape[0]);
}

BOOST_FIXTURE_TEST_CASE(set_fixed_domain_bad_shape, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    space_charge.update_domain(bunch);
    Rectangular_grid_domain_sptr domain_sptr = space_charge.get_domain_sptr();
    std::vector<int > bad_shape(domain_sptr->get_grid_shape());
    // not bad yet
    Rectangular_grid_domain_sptr not_bad_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    space_charge.set_fixed_domain(not_bad_domain_sptr);

    bool caught_error;
    // bad in 0
    caught_error = false;
    ++bad_shape[0];
    Rectangular_grid_domain_sptr bad0_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    try {
        space_charge.set_fixed_domain(bad0_domain_sptr);
    }
    catch (std::runtime_error &) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);

    // bad in 1
    caught_error = false;
    --bad_shape[0];
    ++bad_shape[1];
    Rectangular_grid_domain_sptr bad1_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    try {
        space_charge.set_fixed_domain(bad1_domain_sptr);
    }
    catch (std::runtime_error &) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);

    // bad in 2
    caught_error = false;
    --bad_shape[1];
    ++bad_shape[2];
    Rectangular_grid_domain_sptr bad2_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    try {
        space_charge.set_fixed_domain(bad2_domain_sptr);
    }
    catch (std::runtime_error &) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(get_local_charge_density, Toy_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    std::vector<int > grid_shape_zyx(3);
    grid_shape_zyx[0] = grid_shape[2];
    grid_shape_zyx[1] = grid_shape[1];
    grid_shape_zyx[2] = grid_shape[0];
    Rectangular_grid_domain_sptr domain_sptr = Rectangular_grid_domain_sptr(
            new Rectangular_grid_domain(physical_size, physical_offset,
                    grid_shape_zyx, false));
    space_charge.set_fixed_domain(domain_sptr);
    std::vector<int > center_zyx(3);
    center_zyx[0] = grid_shape_zyx[0] / 2;
    center_zyx[1] = grid_shape_zyx[1] / 2;
    center_zyx[2] = grid_shape_zyx[2] / 2;

    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;
    Rectangular_grid_sptr local_charge_density(
            space_charge.get_local_charge_density(bunch));

    for (int i = center_zyx[0] - 1; i < center_zyx[0] + 1; ++i) {
        for (int j = center_zyx[1] - 1; j < center_zyx[1] + 1; ++j) {
            for (int k = center_zyx[2] - 1; k < center_zyx[2] + 1; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
            }
        }
    }
    multi_array_check_equal(local_charge_density->get_grid_points(), expected,
            100 * tolerance);
}

BOOST_AUTO_TEST_CASE(get_global_charge_density2_reduce_scatter)
{
    Commxx_sptr comm_sptr(new Commxx);

    std::vector<int > grid_shape(3);
    grid_shape[0] = 5;
    grid_shape[1] = 6;
    grid_shape[2] = 16;
    std::vector<int > grid_shape_zyx(3);
    grid_shape_zyx[0] = grid_shape[2];
    grid_shape_zyx[1] = grid_shape[1];
    grid_shape_zyx[2] = grid_shape[0];
    std::vector<double > physical_size(3);
    physical_size[0] = 1.0;
    physical_size[1] = 1.0;
    physical_size[2] = 1.0;
    std::vector<double > physical_offset(3);
    physical_offset[0] = 0.0;
    physical_offset[1] = 0.0;
    physical_offset[2] = 0.0;

    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    Rectangular_grid_domain_sptr domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape_zyx, false));
    space_charge.set_fixed_domain(domain_sptr);
    Rectangular_grid_sptr local_rho(new Rectangular_grid(domain_sptr));
    for (int i = 0; i < grid_shape_zyx[0]; ++i) {
        for (int j = 0; j < grid_shape_zyx[1]; ++j) {
            for (int k = 0; k < grid_shape_zyx[2]; ++k) {
                local_rho->get_grid_points()[i][j][k] = i + 100 * j + 1000 * k;
            }
        }
    }
    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2_reduce_scatter(*local_rho, comm_sptr); // [C/m^3]
    std::vector<int > nondoubled_shape(
            local_rho->get_domain().get_grid_shape());
    std::vector<int > doubled_shape(rho2->get_domain().get_grid_shape());
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[0],
            doubled_shape[0]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[1],
            doubled_shape[1]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[2],
            doubled_shape[2]);

    int lower = rho2->get_lower();
    int upper = rho2->get_upper();
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            for (int k = 0; k < doubled_shape[2]; ++k) {
                if ((i < nondoubled_shape[0]) && (j < nondoubled_shape[1])
                        && (k < nondoubled_shape[2])) {
                    BOOST_CHECK_CLOSE(rho2->get_grid_points()[i][j][k]*
                            rho2->get_normalization(),
                            local_rho->get_grid_points()[i][j][k]*
                            local_rho->get_normalization()*comm_sptr->get_size(),
                            tolerance);
                } else {
                    BOOST_CHECK_EQUAL(rho2->get_grid_points()[i][j][k], 0.0);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(get_global_charge_density2_allreduce)
{
    Commxx_sptr comm_sptr(new Commxx);

    std::vector<int > grid_shape(3);
    grid_shape[0] = 5;
    grid_shape[1] = 6;
    grid_shape[2] = 16;
    std::vector<int > grid_shape_zyx(3);
    grid_shape_zyx[0] = grid_shape[2];
    grid_shape_zyx[1] = grid_shape[1];
    grid_shape_zyx[2] = grid_shape[0];
    std::vector<double > physical_size(3);
    physical_size[0] = 1.0;
    physical_size[1] = 1.0;
    physical_size[2] = 1.0;
    std::vector<double > physical_offset(3);
    physical_offset[0] = 0.0;
    physical_offset[1] = 0.0;
    physical_offset[2] = 0.0;

    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    Rectangular_grid_domain_sptr domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape_zyx, false));
    space_charge.set_fixed_domain(domain_sptr);
    Rectangular_grid_sptr local_rho(new Rectangular_grid(domain_sptr));
    // local_rho will be overwritten by the in-place allreduce, so we
    // keep a copy for comparison.
    Rectangular_grid_sptr local_rho_orig(new Rectangular_grid(domain_sptr));
    for (int i = 0; i < grid_shape_zyx[0]; ++i) {
        for (int j = 0; j < grid_shape_zyx[1]; ++j) {
            for (int k = 0; k < grid_shape_zyx[2]; ++k) {
                local_rho->get_grid_points()[i][j][k] = i + 100 * j + 1000 * k;
                local_rho_orig->get_grid_points()[i][j][k] = i + 100 * j + 1000
                        * k;
            }
        }
    }
    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2_allreduce(*local_rho, comm_sptr); // [C/m^3]
    std::vector<int > nondoubled_shape(
            local_rho->get_domain().get_grid_shape());
    std::vector<int > doubled_shape(rho2->get_domain().get_grid_shape());
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[0],
            doubled_shape[0]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[1],
            doubled_shape[1]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[2],
            doubled_shape[2]);

    int lower = rho2->get_lower();
    int upper = rho2->get_upper();
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            for (int k = 0; k < doubled_shape[2]; ++k) {
                if ((i < nondoubled_shape[0]) && (j < nondoubled_shape[1])
                        && (k < nondoubled_shape[2])) {
                    BOOST_CHECK_CLOSE(rho2->get_grid_points()[i][j][k]*
                            rho2->get_normalization(),
                            local_rho_orig->get_grid_points()[i][j][k]*
                            local_rho_orig->get_normalization()*comm_sptr->get_size(),
                            tolerance);
                } else {
                    BOOST_CHECK_EQUAL(rho2->get_grid_points()[i][j][k], 0.0);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(get_global_charge_density2_simple)
{
    Commxx_sptr comm_sptr(new Commxx);

    std::vector<int > grid_shape(3);
    grid_shape[0] = 5;
    grid_shape[1] = 6;
    grid_shape[2] = 16;
    std::vector<int > grid_shape_zyx(3);
    grid_shape_zyx[0] = grid_shape[2];
    grid_shape_zyx[1] = grid_shape[1];
    grid_shape_zyx[2] = grid_shape[0];
    std::vector<double > physical_size(3);
    physical_size[0] = 1.0;
    physical_size[1] = 1.0;
    physical_size[2] = 1.0;
    std::vector<double > physical_offset(3);
    physical_offset[0] = 0.0;
    physical_offset[1] = 0.0;
    physical_offset[2] = 0.0;

    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    Rectangular_grid_domain_sptr domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape_zyx, false));
    space_charge.set_fixed_domain(domain_sptr);
    Rectangular_grid_sptr local_rho(new Rectangular_grid(domain_sptr));
    // local_rho will be overwritten by the in-place allreduce, so we
    // keep a copy for comparison.
    Rectangular_grid_sptr local_rho_orig(new Rectangular_grid(domain_sptr));
    for (int i = 0; i < grid_shape_zyx[0]; ++i) {
        for (int j = 0; j < grid_shape_zyx[1]; ++j) {
            for (int k = 0; k < grid_shape_zyx[2]; ++k) {
                local_rho->get_grid_points()[i][j][k] = i + 100 * j + 1000 * k;
                local_rho_orig->get_grid_points()[i][j][k] = i + 100 * j + 1000
                        * k;
            }
        }
    }
    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2(*local_rho, comm_sptr); // [C/m^3]
    std::vector<int > nondoubled_shape(
            local_rho->get_domain().get_grid_shape());
    std::vector<int > doubled_shape(rho2->get_domain().get_grid_shape());
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[0],
            doubled_shape[0]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[1],
            doubled_shape[1]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain().get_grid_shape()[2],
            doubled_shape[2]);

    int lower = rho2->get_lower();
    int upper = rho2->get_upper();
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            for (int k = 0; k < doubled_shape[2]; ++k) {
                if ((i < nondoubled_shape[0]) && (j < nondoubled_shape[1])
                        && (k < nondoubled_shape[2])) {
                    BOOST_CHECK_CLOSE(rho2->get_grid_points()[i][j][k]*
                            rho2->get_normalization(),
                            local_rho_orig->get_grid_points()[i][j][k]*
                            local_rho_orig->get_normalization()*comm_sptr->get_size(),
                            tolerance);
                } else {
                    BOOST_CHECK_EQUAL(rho2->get_grid_points()[i][j][k], 0.0);
                }
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_green_fn2_pointlike, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    space_charge.update_domain(bunch);
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike());
    MArray3d_ref G2_a(G2->get_grid_points());
    double norm = G2->get_normalization();
    int jmirror, kmirror;
    double dz, dy, dx;
    const double coeff = 2.8;
    double G000 = coeff / std::min(G2->get_domain().get_cell_size()[0],
            std::min(G2->get_domain().get_cell_size()[1],
                    G2->get_domain().get_cell_size()[2]));

    for (int i = G2->get_lower(); i < G2->get_upper(); ++i) {
        int i_dz;
        if (i < G2->get_domain().get_grid_shape()[0] / 2) {
            i_dz = i;
        } else {
            i_dz = G2->get_domain().get_grid_shape()[0] - i;
        }
        dz = i_dz * G2->get_domain().get_cell_size()[0];
        for (int j = 0; j < G2->get_domain().get_grid_shape()[1] / 2; ++j) {
            dy = j * G2->get_domain().get_cell_size()[1];
            jmirror = G2->get_domain().get_grid_shape()[1] - j;
            if (jmirror == G2->get_domain().get_grid_shape()[1]) {
                jmirror = j;
            }
            for (int k = 0; k < G2->get_domain().get_grid_shape()[2] / 2; ++k) {
                dx = k * G2->get_domain().get_cell_size()[2];
                kmirror = G2->get_domain().get_grid_shape()[2] - k;
                if (kmirror == G2->get_domain().get_grid_shape()[2]) {
                    kmirror = k;
                }
                double G;
                if ((i == 0) && (j == 0) && (k == 0)) {
                    G = G000;
                } else {
                    G = 1 / std::sqrt(dx * dx + dy * dy + dz * dz);
                }
                BOOST_CHECK_CLOSE(G2_a[i][j][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][jmirror][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][jmirror][kmirror]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][j][kmirror]*norm, G, tolerance);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_green_fn2_pointlike_periodic, Ellipsoidal_bunch_fixture)
{
    double z_period = 10.0 * stdz;
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, true, true,
            z_period);
    space_charge.update_domain(bunch);
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike());
    MArray3d_ref G2_a(G2->get_grid_points());
    double norm = G2->get_normalization();
    int jmirror, kmirror;
    double dz, dy, dx;
    const double coeff = 2.8;
    double G000 = coeff / std::min(G2->get_domain().get_cell_size()[0],
            std::min(G2->get_domain().get_cell_size()[1],
                    G2->get_domain().get_cell_size()[2]));
    const int num_images = 8;
    for (int i = G2->get_lower(); i < G2->get_upper(); ++i) {
        int i_dz;
        if (i < G2->get_domain().get_grid_shape()[0] / 2) {
            i_dz = i;
        } else {
            i_dz = G2->get_domain().get_grid_shape()[0] - i;
        }
        dz = i_dz * G2->get_domain().get_cell_size()[0];
        for (int j = 0; j < G2->get_domain().get_grid_shape()[1] / 2; ++j) {
            dy = j * G2->get_domain().get_cell_size()[1];
            jmirror = G2->get_domain().get_grid_shape()[1] - j;
            if (jmirror == G2->get_domain().get_grid_shape()[1]) {
                jmirror = j;
            }
            for (int k = 0; k < G2->get_domain().get_grid_shape()[2] / 2; ++k) {
                dx = k * G2->get_domain().get_cell_size()[2];
                kmirror = G2->get_domain().get_grid_shape()[2] - k;
                if (kmirror == G2->get_domain().get_grid_shape()[2]) {
                    kmirror = k;
                }
                double G;
                if ((i == 0) && (j == 0) && (k == 0)) {
                    G = G000;
                } else {
                    G = 1 / std::sqrt(dx * dx + dy * dy + dz * dz);
                }
                for (int image = -num_images; image <= num_images; ++image) {
                    double dz_image = dz + image * z_period;
                    if (image != 0) {
                        if ((j == 0) && (k == 0) && (std::abs(dz_image)
                                < 1.0e-9)) {
                            G += G000;
                        } else {
                            G += 1.0 / std::sqrt(dx * dx + dy * dy + dz_image
                                    * dz_image);
                        }
                    }
                }
                BOOST_CHECK_CLOSE(G2_a[i][j][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][jmirror][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][jmirror][kmirror]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][j][kmirror]*norm, G, tolerance);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_green_fn2_no_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape);
    bool caught_error = false;
    try {
        Distributed_rectangular_grid_sptr G2(
                space_charge.get_green_fn2_pointlike());
    }
    catch (std::runtime_error &) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);
}

//// solver tests in test_space_charge_3d_open_hockney(2-4).cc
//
//void
//simple_populate(Bunch & bunch, Random_distribution & distribution)
//{
//    MArray2d covariances(boost::extents[6][6]);
//    MArray1d means(boost::extents[6]);
//    for (int i = 0; i < 6; ++i) {
//        means[i] = 0.0;
//        for (int j = i; j < 6; ++j) {
//            covariances[i][j] = 0.0;
//        }
//    }
//    // This bunch shape is contrived to make longitudinal kicks be comparable
//    // to transverse kicks
//    double stdx = 1.1e-3;
//    double stdy = 2.3e-3;
//    double stdz = 3.5e-7;
//    covariances[0][0] = stdx * stdx;
//    covariances[2][2] = stdy * stdy;
//    covariances[4][4] = stdz * stdz;
//    covariances[1][1] = covariances[3][3] = covariances[5][5] = 1.0e-3;
//    populate_6d(distribution, bunch, means, covariances);
//
//}
//
//BOOST_FIXTURE_TEST_CASE(apply_full, Ellipsoidal_bunch_fixture)
//{
//    simple_populate(bunch, distribution);
//    Bunch original_bunch(bunch);
//    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, true);
//    const double time_fraction = 1.0;
//    Step dummy_step(time_fraction);
//    const double time_step = 0.3;
//    space_charge.apply(bunch, time_step, dummy_step);
//
//    double total_x_kick2 = 0.0;
//    double total_y_kick2 = 0.0;
//    double total_p_kick2 = 0.0;
//    for (int i = 0; i < bunch.get_local_num(); ++i) {
//        double kick;
//        kick = bunch.get_local_particles()[i][Bunch::xp]
//                - original_bunch.get_local_particles()[i][Bunch::xp];
//        total_x_kick2 += kick * kick;
//        kick += bunch.get_local_particles()[i][Bunch::yp]
//                - original_bunch.get_local_particles()[i][Bunch::yp];
//        total_y_kick2 += kick * kick;
//        kick = bunch.get_local_particles()[i][Bunch::dpop]
//                - original_bunch.get_local_particles()[i][Bunch::dpop];
//        total_p_kick2 += kick * kick;
//    }
//    double avg_x_kick2 = total_x_kick2 / bunch.get_local_num();
//    double avg_y_kick2 = total_y_kick2 / bunch.get_local_num();
//    double avg_p_kick2 = total_p_kick2 / bunch.get_local_num();
//
//    const double rough_tolerance = 50.0;
//    BOOST_CHECK_CLOSE(avg_x_kick2, 2.6e5, rough_tolerance);
//    BOOST_CHECK_CLOSE(avg_y_kick2, 3.4e5, rough_tolerance);
//    BOOST_CHECK_CLOSE(avg_p_kick2, 7.0e7, rough_tolerance);
//}
//
//BOOST_FIXTURE_TEST_CASE(apply_transverse, Ellipsoidal_bunch_fixture)
//{
//    simple_populate(bunch, distribution);
//    Bunch original_bunch(bunch);
//    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, false);
//    const double time_fraction = 1.0;
//    Step dummy_step(time_fraction);
//    const double time_step = 0.3;
//    space_charge.apply(bunch, time_step, dummy_step);
//
//    double total_x_kick2 = 0.0;
//    double total_y_kick2 = 0.0;
//    double total_p_kick2 = 0.0;
//    for (int i = 0; i < bunch.get_local_num(); ++i) {
//        double kick;
//        kick = bunch.get_local_particles()[i][Bunch::xp]
//                - original_bunch.get_local_particles()[i][Bunch::xp];
//        total_x_kick2 += kick * kick;
//        kick += bunch.get_local_particles()[i][Bunch::yp]
//                - original_bunch.get_local_particles()[i][Bunch::yp];
//        total_y_kick2 += kick * kick;
//        kick = bunch.get_local_particles()[i][Bunch::dpop]
//                - original_bunch.get_local_particles()[i][Bunch::dpop];
//        total_p_kick2 += kick * kick;
//    }
//    double avg_x_kick2 = total_x_kick2 / bunch.get_local_num();
//    double avg_y_kick2 = total_y_kick2 / bunch.get_local_num();
//    double avg_p_kick2 = total_p_kick2 / bunch.get_local_num();
//
//    const double rough_tolerance = 50.0;
//    BOOST_CHECK_CLOSE(avg_x_kick2, 2.6e5, rough_tolerance);
//    BOOST_CHECK_CLOSE(avg_y_kick2, 3.4e5, rough_tolerance);
//    BOOST_CHECK_CLOSE(avg_p_kick2, 7.1e4, rough_tolerance);
//}
#endif
