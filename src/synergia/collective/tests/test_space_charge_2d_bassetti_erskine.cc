#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/space_charge_2d_bassetti_erskine.h"
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
#include "synergia/utils/hdf5_writer.h"
#include "synergia/bunch/core_diagnostics.h"
#include "gaussian_charge_density.h"
//#include "space_charge_bunch_fixtures.h"
//#include "synergia/utils/simple_timer.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.7e11;
const int total_num = 10000;
const double total_energy = 125.0;

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Space_charge_2d_bassetti_erskine space_charge;
    double sigma_x = 3.5e-3;
    double sigma_y = 2.3e-3;
    double sigma_cdt = 1.2e-2;
    space_charge.set_sigma(sigma_x, sigma_y, sigma_cdt);
}

void
simple_populate(Bunch & bunch, Random_distribution & distribution)
{
    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] = 0.0;
        for (int j = i; j < 6; ++j) {
            covariances[i][j] = 0.0;
        }
    }
    // This bunch shape is contrived to make longitudinal kicks be comparable
    // to transverse kicks
    double stdx = 1.1e-3;
    double stdy = 2.3e-3;
    double stdz = 3.5e-7;
    covariances[0][0] = stdx * stdx;
    covariances[2][2] = stdy * stdy;
    covariances[4][4] = stdz * stdz;
    covariances[1][1] = covariances[3][3] = covariances[5][5] = 1.0e-3;
    populate_6d(distribution, bunch, means, covariances);

}


//const double stdx = 0.0002;  // for round rod
const double stdx = 0.000001;  // for round rod
const double stdy = 0.0001;  // for elliptical rod
const int rod_num_particles = 1 + 5001*8; // 1 + (odd number) * 8 + 2 to expand the domein on both sides
const double rod_lowgamma = 61.0/60.0;
const double rod_highgamma = 61.0/11.0;
const double rod_real_num = 5.0e9;
const double rod_length = 0.1;
// the RMS of the distribution is rod_radius/sqrt(2)
//const double rod_radius = 1.0e-3; //radius of rod particles
const double rod_probe = 0.001;  // radius of the probe particle

struct Round_rod_bunch_fixture_lowgamma
{
    Round_rod_bunch_fixture_lowgamma() :
        four_momentum(mass, mass*rod_lowgamma), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                rod_num_particles, rod_real_num,
                            comm_sptr, rod_length)
    {
        BOOST_TEST_MESSAGE("setup Round_rod bunch fixture lowgamma");
        const double rod_radius = stdx * std::sqrt(2.0);
        bunch.set_sort_period(-1);
        MArray2d_ref local_particles(bunch.get_local_particles());
        // a ring of 8 particles around each longitudinal location
        int num_longitudinal = (rod_num_particles-1)/8;
        double dz = rod_length/(num_longitudinal-1);

        double r2o2 = std::sqrt(2.0)/2.0;

        double z = -rod_length/2.0;
        for (int i=1; i<rod_num_particles; i+=8, z+=dz) {
            local_particles[i][Bunch::x] = rod_radius;
            local_particles[i][Bunch::y] = 0.0;

            local_particles[i+1][Bunch::x] = rod_radius*r2o2;
            local_particles[i+1][Bunch::y] = rod_radius*r2o2;

            local_particles[i+2][Bunch::x] = 0.0;
            local_particles[i+2][Bunch::y] = rod_radius;

            local_particles[i+3][Bunch::x] = -rod_radius*r2o2;
            local_particles[i+3][Bunch::y] =  rod_radius*r2o2;

            local_particles[i+4][Bunch::x] = -rod_radius;
            local_particles[i+4][Bunch::y] = 0.0;

            local_particles[i+5][Bunch::x] = -rod_radius*r2o2;
            local_particles[i+5][Bunch::y] = -rod_radius*r2o2;

            local_particles[i+6][Bunch::x] = 0.0;
            local_particles[i+6][Bunch::y] = -rod_radius;

            local_particles[i+7][Bunch::x] = rod_radius*r2o2;
            local_particles[i+7][Bunch::y] = -rod_radius*r2o2;

            double rod_beta = reference_particle.get_beta();
            for (int j=i; j<i+8; ++j) {
                local_particles[j][Bunch::cdt] = z/rod_beta;
                local_particles[j][Bunch::xp] = 0.0;
                local_particles[j][Bunch::yp] = 0.0;
                local_particles[j][Bunch::dpop] = 0.0;
                local_particles[j][Bunch::id] = j;
            }
        }

        // Bassetti-Erskine solver uses calculated RMS so particle can be anywhere
        local_particles[0][Bunch::x] = rod_probe;
        local_particles[0][Bunch::y] = 0.0;
        local_particles[0][Bunch::xp] = 0.0;
        local_particles[0][Bunch::yp] = 0.0;
        local_particles[0][Bunch::cdt] = 0.0;
        local_particles[0][Bunch::dpop] = 0.0;
        local_particles[0][Bunch::id] = 0.0;

    }

    ~Round_rod_bunch_fixture_lowgamma()
    {
        BOOST_TEST_MESSAGE("tear down Rod bunch fixture lowgamma");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
};

BOOST_FIXTURE_TEST_CASE(apply_bassetti_erskine_round, Round_rod_bunch_fixture_lowgamma)
{
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    Space_charge_2d_bassetti_erskine space_charge;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    logger << "probe radius: " << probe_radius << std::endl;

    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[1] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[1] << std::endl;
    double stdcdt = stds[2];

    space_charge.apply(bunch, time_step, dummy_step, verbosity, logger);
    bunch.convert_to_state(Bunch::fixed_z_lab);

    // Rod of charge Q over length L of gaussian with RMS sigma
    // Fraction of charge at radius R is 1 - exp(R/sigma)
    // E field at radius r $$ E = \frac{1}{2 \pi \epsilon_0} \frac{Q}{L} \frac{1}{r} $$
    // B field at radius r $$ E = \frac{\mu_0}{2 \pi } \frac{Q v}{L} \frac{1}{r} $$
    // Net EM force on electric+magnetic on probe of charge q from E-B cancellation
    // $$ F = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{1}{r}
    // travel over distance D at velocity v
    // \frac{\Delta p}{p} = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{D}{m v^2} \frac{1}{r}
    // convert to usual units
    // \frac{\Delta p}{p} = \frac{2 N r_p}{L \beta^2 \gamma^3} \frac{D}{r}

    double L = bunch.get_z_period_length();
    double fractional_charge = (1.0 - std::exp(-probe_radius*probe_radius/(2.0*stdx*stdx))/(std::sqrt(2.0*mconstants::pi)*stdcdt));
    double N = bunch.get_real_num()*fractional_charge;

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(L*betagamma*betagamma*gamma)) *
            (step_length/bunch.get_local_particles()[0][Bunch::x]);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);

}

#if 0
BOOST_FIXTURE_TEST_CASE(efield_particles, Spherical_bunch_fixture)
{
    double sigin[3];
    sigin[0] = sigmax;
    sigin[1] = sigmay;
    sigin[2] = sigmaz;
    Space_charge_2d_bassetti_erskine space_charge_bs;
    space_charge_bs.set_sigma(sigin);
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
    * pconstants::e;

    for (int component = 0; component < 1; ++component) {
        double max_fractional_error = -2.0;
        double min_fractional_error = 2.0;
        int npoints = 100;
        double x, y, z;
        double xmin = -0.005, xmax = 0.005;
        double ymin = -0.005, ymax = 0.005;
        double zmin = -0.005, zmax = 0.005;
        for (int i = 0; i <= npoints; ++i) {
            x = xmin + i * (xmax - xmin) / npoints;
//            for (int j = 0; j <= npoints; ++j) {
//                y = ymin + j * (ymax - ymin) / npoints;
//                for (int k = 0; k <= npoints; ++k) {
//                    z = zmin + k * (zmax - zmin) / npoints;
            y = 0, z = 0;

            double r = std::sqrt(x * x + y * y + z * z);
            double var = x;
            double efield_exact, efield_calc;

            double line_charge_density = Q * exp(-z * z /(2.0
                            * sigin[2] * sigin[2])) / (sqrt(2.0
                            * mconstants::pi) * sigin[2]);
            if (component == 0) {
                var = x;
                efield_exact = gaussian_electric_field_component2(Q,
                        x, y, z, sigin[0], sigin[1], sigin[2],
                        var);
                efield_calc = space_charge_bs.normalized_efield(x,y)[0]
                / (2.0 * mconstants::pi * pconstants::epsilon0)
                * line_charge_density;
            } else if (component == 1) {
                var = y;
                efield_exact = gaussian_electric_field_component2(Q,
                        x, y, z, sigin[0], sigin[1], sigin[2],
                        var);
                efield_calc = space_charge_bs.normalized_efield(x,y)[1]
                / (2.0 * mconstants::pi * pconstants::epsilon0)
                * line_charge_density;
            }

            const double tiny = 1.0e-8;

            double fractional_error = (efield_calc - efield_exact)
            / efield_exact;
            if (fractional_error > max_fractional_error) {
                max_fractional_error = fractional_error;
            }
            if (fractional_error < min_fractional_error) {
                min_fractional_error = fractional_error;
            }
#if 1
            std::cout << x << "  " << y << "  " << z << "  "
            << efield_exact << "  "
            << efield_calc << "  "
            << fractional_error << std::endl;
#endif
//                }
//            }
        }
        std::cout << "max_fractional_error = " << max_fractional_error
        << std::endl;
        std::cout << "min_fractional_error = " << min_fractional_error
        << std::endl;
        const double field_tolerance[] = {2.0, 2.0};
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}
#endif

//BOOST_FIXTURE_TEST_CASE(efield_particles, Spherical_bunch_fixture_2d)
//{
//    std::vector<int > grid_shape_xyz(3);
//    grid_shape_xyz[0] = 128;
//    grid_shape_xyz[1] = 128;
//    grid_shape_xyz[2] = 64;
//
//    double t;
//    t = simple_timer_current();
//    // Bassetti-Erskine solver
//    Space_charge_2d_bassetti_erskine space_charge_bs;
//    MArray1d mean(Diagnostics::calculate_mean(bunch));
//    MArray1d std(Diagnostics::calculate_std(bunch, mean));
//    double sigin[3];
//    sigin[0] = std[Bunch::x];
//    sigin[1] = std[Bunch::y];
//    sigin[2] = std[Bunch::z];
//    space_charge_bs.set_sigma(sigin);
//    t = simple_timer_show(t, "bs construct");
//    // 2.5D Open Hockney solver
//    Space_charge_2d_open_hockney space_charge_2d(comm, grid_shape_xyz);
//    //double t;
//    //t = simple_timer_current();
//    Rectangular_grid_sptr local_rho(
//            space_charge_2d.get_local_charge_density(bunch)); // [C/m^3]
//    Distributed_rectangular_grid_sptr rho2 =
//            space_charge_2d.get_global_charge_density2(*local_rho); // [C/m^3]
//    //t = simple_timer_show(t, "charge_deposition");
//    Distributed_rectangular_grid_sptr
//            G2(space_charge_2d.get_green_fn2_pointlike()); // [1/m]
//    Distributed_rectangular_grid_sptr local_force2(
//            space_charge_2d.get_local_force2(*rho2, *G2)); // [N]
//    t = simple_timer_show(t, "2.5D construct");
//    // 3D Open Hockney solver
//    Space_charge_3d_open_hockney space_charge_3d(comm, grid_shape_xyz);
//    Rectangular_grid_sptr local_rho_3d(
//            space_charge_3d.get_local_charge_density(bunch)); // [C/m^3]
//    Distributed_rectangular_grid_sptr rho2_3d(
//            space_charge_3d.get_global_charge_density2(*local_rho_3d)); // [C/m^3]
//    Distributed_rectangular_grid_sptr
//            G2_3d(space_charge_3d.get_green_fn2_linear()); // [1/m]
//    Distributed_rectangular_grid_sptr phi2(space_charge_3d.get_scalar_field2(
//            *rho2_3d, *G2_3d)); // [V]
//    Distributed_rectangular_grid_sptr phi(space_charge_3d.extract_scalar_field(
//            *phi2));
//    phi->fill_guards();
//    t = simple_timer_show(t, "3D construct");
//
//    double Q = bunch.get_real_num() * bunch.get_particle_charge()
//            * pconstants::e;
//    double q = bunch.get_particle_charge() * pconstants::e;
//    std::vector<int > nondoubled_shape(
//             space_charge_2d.get_domain().get_grid_shape());
//    std::vector<int > doubled_shape(
//             space_charge_2d.get_doubled_domain_sptr()->get_grid_shape());
//    double hz = space_charge_2d.get_domain().get_cell_size()[2];
//
//    for (int component = 0; component < 1; ++component) {
//        double max_fractional_error = -2.0;
//        double min_fractional_error = 2.0;
//        Distributed_rectangular_grid_sptr local_En(
//                space_charge_3d.get_electric_field_component(*phi, component));
//        Rectangular_grid_sptr En(
//                space_charge_3d.get_global_electric_field_component(*local_En));
//        for (int i = nondoubled_shape[0] / 2; i < 3 * nondoubled_shape[0] / 2;
//                ++i) {
////            for (int j = nondoubled_shape[1] / 2; j < 3 * nondoubled_shape[1]
////                    / 2; ++j) {
////                for (int k = 0; k < doubled_shape[2]; ++k) {
////                    int i = nondoubled_shape[0];
//                    int j = nondoubled_shape[1];
//                    int k = nondoubled_shape[2] / 2;
//                    double x, y, z;
//                    space_charge_2d.get_doubled_domain_sptr()->get_cell_coordinates(i, j, k, x, y, z);
//
//                    int ip = i - nondoubled_shape[0] / 2;
//                    int jp = nondoubled_shape[1] / 2;
//                    int kp = nondoubled_shape[2] / 2;
//                    double xp, yp , zp;
//                    space_charge_3d.get_domain().get_cell_coordinates(
//                            kp, jp, ip, zp, yp, xp);
//
//                    double r = std::sqrt(x * x + y * y + z * z);
//                    double var = x;
//                    double efield_exact, efield_calc3, efield_calc4;
//                    double efield_calc5, efield_calc6;
//
//                    double line_charge_density1 = Q * exp(-z * z /(2.0
//                            * sigin[2] * sigin[2])) / (sqrt(2.0
//                            * mconstants::pi) * sigin[2]);
//                    double line_charge_density2 =  bunch.get_real_num()
//                            / bunch.get_total_num()
//                            * bunch.get_particle_charge() * pconstants::e
//                            * rho2->get_grid_points_1d()[k] / hz;
//
//                    efield_calc3 = space_charge_bs.normalized_efield(x,y)[0]
//                            / (2.0 * mconstants::pi * pconstants::epsilon0)
//                            * line_charge_density1;
//                    efield_calc4 = space_charge_bs.normalized_efield(x,y)[0]
//                            / (2.0 * mconstants::pi * pconstants::epsilon0)
//                            * line_charge_density2;
//                    efield_calc5 = rho2->get_grid_points_1d()[k]
//                                * local_force2->get_grid_points_2dc()[i][j].real()
//                                * local_force2->get_normalization() / q;
//                    efield_calc6 = En->get_grid_points()[kp][jp][ip]
//                                * En->get_normalization();
//
//                    efield_exact
//                            = gaussian_electric_field_component2(Q, x, y, z,
//                                    sigin[0], sigin[1], sigin[2], var);
//
//                    const double tiny = 1.0e-8;
//
//                    double fractional_error1 = (efield_calc3 - efield_exact)
//                                / efield_exact;
//                    double fractional_error2 = (efield_calc4 - efield_exact)
//                                / efield_exact;
//                    double fractional_error3 = (efield_calc5 - efield_exact)
//                                / efield_exact;
//                    double fractional_error4 = (efield_calc6 - efield_exact)
//                                / efield_exact;
//                    if (fractional_error1 > max_fractional_error) {
//                        max_fractional_error = fractional_error1;
//                    }
//                    if (fractional_error1 < min_fractional_error) {
//                        min_fractional_error = fractional_error1;
//                    }
//
//                    #if 1
//                    std::cout << x << "  " << y << "  " << z << "  "
//                            << efield_exact << "  "
//                            << efield_calc3 << "  "
//                            << efield_calc4 << "  "
//                            << efield_calc5 << "  "
//                            << efield_calc6 << "  "
//                            << fractional_error1 << "  "
//                            << fractional_error2 << "  "
//                            << fractional_error3 << "  "
//                            << fractional_error4 << std::endl;
//                    #endif
////                }
////            }
//        }
//        std::cout << "max_fractional_error = " << max_fractional_error
//                << std::endl;
//        std::cout << "min_fractional_error = " << min_fractional_error
//                << std::endl;
//        const double field_tolerance[] = { 2.0, 2.0 };
//        //BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
//        //BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
//    }
//}
