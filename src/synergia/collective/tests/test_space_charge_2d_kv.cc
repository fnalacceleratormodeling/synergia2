#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/space_charge_2d_kv.h"
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
    Space_charge_2d_kv space_charge;
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


const double stdx = 0.0005;
const double stdy = 0.0002;  // for elliptical rod
const int rod_num_particles = 4 + 100001*8; // 4 + (odd number) * 8 to add 4 test particles that don't
                                          // change mean or std
const double rod_lowgamma = 61.0/60.0;
const double rod_highgamma = 61.0/11.0;
const double rod_real_num = 5.0e9;
const double rod_length = 0.1;
// the RMS of the distribution is rod_radius/sqrt(2)
const double rod_probe = 10.0*stdx;  // far outside of rod
struct Elliptical_rod_bunch_fixture_lowgamma
{
    Elliptical_rod_bunch_fixture_lowgamma() :
        four_momentum(mass, mass*rod_lowgamma), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                rod_num_particles, rod_real_num,
                            comm_sptr, rod_length)
    {
        BOOST_TEST_MESSAGE("setup Round_rod bunch fixture lowgamma");
        const double rod_xradius = stdx * std::sqrt(2.0);
        const double rod_yradius = stdy * std::sqrt(2.0);
        bunch.set_sort_period(-1);
        MArray2d_ref local_particles(bunch.get_local_particles());
        // a ring of 8 particles around each longitudinal location
        int num_longitudinal = (rod_num_particles-1)/8;
        double dz = rod_length/(num_longitudinal-1);

        double r2o2 = std::sqrt(2.0)/2.0;

        double z = -rod_length/2.0;
        for (int i=4; i<rod_num_particles; i+=8, z+=dz) {
            local_particles[i][Bunch::x] = rod_xradius;
            local_particles[i][Bunch::y] = 0.0;

            local_particles[i+1][Bunch::x] = rod_xradius*r2o2;
            local_particles[i+1][Bunch::y] = rod_yradius*r2o2;

            local_particles[i+2][Bunch::x] = 0.0;
            local_particles[i+2][Bunch::y] = rod_yradius;

            local_particles[i+3][Bunch::x] = -rod_xradius*r2o2;
            local_particles[i+3][Bunch::y] =  rod_yradius*r2o2;

            local_particles[i+4][Bunch::x] = -rod_xradius;
            local_particles[i+4][Bunch::y] = 0.0;

            local_particles[i+5][Bunch::x] = -rod_xradius*r2o2;
            local_particles[i+5][Bunch::y] = -rod_yradius*r2o2;

            local_particles[i+6][Bunch::x] = 0.0;
            local_particles[i+6][Bunch::y] = -rod_yradius;

            local_particles[i+7][Bunch::x] = rod_xradius*r2o2;
            local_particles[i+7][Bunch::y] = -rod_yradius*r2o2;

            double rod_beta = reference_particle.get_beta();
            for (int j=i; j<i+8; ++j) {
                local_particles[j][Bunch::cdt] = z/rod_beta;
                local_particles[j][Bunch::xp] = 0.0;
                local_particles[j][Bunch::yp] = 0.0;
                local_particles[j][Bunch::dpop] = 0.0;
                local_particles[j][Bunch::id] = j;
            }
        }

        // Add test particles at += stdx so mean and std are
        // not shifted
        local_particles[0][Bunch::x] = rod_probe;
        local_particles[0][Bunch::y] = 0.0;
        local_particles[0][Bunch::xp] = 0.0;
        local_particles[0][Bunch::yp] = 0.0;
        local_particles[0][Bunch::cdt] = 0.0;
        local_particles[0][Bunch::dpop] = 0.0;
        local_particles[0][Bunch::id] = 0.0;

        local_particles[1][Bunch::x] = -rod_probe;
        local_particles[1][Bunch::y] = 0.0;
        local_particles[1][Bunch::xp] = 0.0;
        local_particles[1][Bunch::yp] = 0.0;
        local_particles[1][Bunch::cdt] = 0.0;
        local_particles[1][Bunch::dpop] = 0.0;
        local_particles[1][Bunch::id] = 0.0;

        local_particles[2][Bunch::x] = 0.0;
        local_particles[2][Bunch::y] = rod_probe;
        local_particles[2][Bunch::xp] = 0.0;
        local_particles[2][Bunch::yp] = 0.0;
        local_particles[2][Bunch::cdt] = 0.0;
        local_particles[2][Bunch::dpop] = 0.0;
        local_particles[2][Bunch::id] = 0.0;

        local_particles[3][Bunch::x] = 0.0;
        local_particles[3][Bunch::y] = -rod_probe;
        local_particles[3][Bunch::xp] = 0.0;
        local_particles[3][Bunch::yp] = 0.0;
        local_particles[3][Bunch::cdt] = 0.0;
        local_particles[3][Bunch::dpop] = 0.0;
        local_particles[3][Bunch::id] = 0.0;
    }

    ~Elliptical_rod_bunch_fixture_lowgamma()
    {
        BOOST_TEST_MESSAGE("tear down Rod bunch fixture lowgamma");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
};

BOOST_FIXTURE_TEST_CASE(apply_kv_elliptical_lowgamma, Elliptical_rod_bunch_fixture_lowgamma)
{
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    bunch.convert_to_state(Bunch::fixed_z_lab);
    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[2] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[2] << std::endl;
    double stdcdt = stds[4];
    logger << "Bunch stdcdt: " << stdcdt << std::endl;

    Space_charge_2d_kv space_charge;
    std::cout << "KV longitudinal: " << space_charge.get_longitudinal() << std::endl;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    logger << "probe radius: " << probe_radius << std::endl;


    space_charge.apply(bunch, time_step, dummy_step, verbosity, logger);
    bunch.convert_to_state(Bunch::fixed_z_lab);

    // Rod of charge Q over length L
    // E field at radius r $$ E = \frac{1}{2 \pi \epsilon_0} \frac{Q}{L}
    // \frac{1}{r} $$
    // B field at radius r $$ E = \frac{\mu_0}{2 \pi } \frac{Q v}{L} \frac{1}{r}
    // $$
    // Net EM force on electric+magnetic on probe of charge q from E-B
    // cancellation
    // $$ F = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{1}{r}
    // travel over distance D at velocity v
    // \frac{\Delta p}{p} = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L}
    // \frac{D}{m v^2} \frac{1}{r}
    // convert to usual units
    // \frac{\Delta p}{p} = \frac{2 N r_p}{L \beta^2 \gamma^3} \frac{D}{r}
    double L = bunch.get_z_period_length();
    double N = bunch.get_real_num();
    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop =
        ((2.0 * N * pconstants::rp) / (L * betagamma * betagamma * gamma)) *
        (step_length / bunch.get_local_particles()[0][Bunch::x]);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[0][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[1][Bunch::xp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[1][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[2][Bunch::yp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[2][Bunch::xp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[3][Bunch::yp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[3][Bunch::xp]), 1.0e-10);
}

#if 0
BOOST_FIXTURE_TEST_CASE(apply_kv_round_gaussian_lowgamma, Round_rod_bunch_fixture_lowgamma)
{
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    bunch.convert_to_state(Bunch::fixed_z_lab);
    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[2] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[2] << std::endl;
    double stdcdt = stds[4];
    logger << "Bunch stdcdt: " << stdcdt << std::endl;

    Space_charge_2d_kv space_charge;
    // set to use gaussian charge density
    space_charge.set_longitudinal(Space_charge_2d_kv::longitudinal_gaussian);
    std::cout << "bassetti-erskine longitudinal: " << space_charge.get_longitudinal() << std::endl;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    logger << "probe radius: " << probe_radius << std::endl;


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
    // volume factor is the charge contained in a cylinder of radius smaller than the probe
    double volume_factor = 1.0 - std::exp(-probe_radius*probe_radius/(2.0*stdx*stdx));
    logger << "volume_factor: " << volume_factor << std::endl;
    // if fractional_charge uses longitudinal normal distribution
    double length_normalization = 1.0/(std::sqrt(2.0*mconstants::pi)*stdcdt*beta);

    logger << "length normalization factor: " << length_normalization << std::endl;
    double N = bunch.get_real_num()*volume_factor * length_normalization;

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(betagamma*betagamma*gamma)) *
            (step_length/bunch.get_local_particles()[0][Bunch::x]);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[0][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[1][Bunch::xp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[1][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[2][Bunch::yp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[2][Bunch::xp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[3][Bunch::yp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[3][Bunch::xp]), 1.0e-10);
}

BOOST_FIXTURE_TEST_CASE(apply_kv_round_offset_lowgamma, Round_rod_bunch_fixture_lowgamma)
{
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    bunch.convert_to_state(Bunch::fixed_z_lab);
    // add offset to bunch
    for (int i=0; i<bunch.get_local_num(); ++i) {
        bunch.get_local_particles()[i][Bunch::x] += 1.0;
        bunch.get_local_particles()[i][Bunch::y] += -2.0;
    }
    
    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[2] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[2] << std::endl;
    double stdcdt = stds[4];
    logger << "Bunch stdcdt: " << stdcdt << std::endl;

    Space_charge_2d_kv space_charge;
    std::cout << "bassetti-erskine longitudinal: " << space_charge.get_longitudinal() << std::endl;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    double probe_offset = bunch.get_local_particles()[0][Bunch::x]-mean[0];
    logger << "probe radius: " << probe_radius << std::endl;
    logger << "probe_offset: " << probe_offset << std::endl;

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
    // volume factor is the charge contained in a cylinder of radius smaller than the probe
    double volume_factor = 1.0 - std::exp(-probe_offset*probe_offset/(2.0*stdx*stdx));
    logger << "volume_factor: " << volume_factor << std::endl;
    // if fractional_charge includes the assumed longitudinal normal distribution
    // double length_normalization = 1.0/(std::sqrt(2.0*mconstants::pi)*stdcdt*beta);
    double length_normalization = 1.0/bunchlen;
    logger << "length normalization factor: " << length_normalization << std::endl;
    double N = bunch.get_real_num()*volume_factor * length_normalization;

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(betagamma*betagamma*gamma)) *
            (step_length/probe_offset);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[0][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[1][Bunch::xp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[1][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[2][Bunch::yp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[2][Bunch::xp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[3][Bunch::yp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[3][Bunch::xp]), 1.0e-10);
}

struct Round_rod_bunch_fixture_highgamma
{
    Round_rod_bunch_fixture_highgamma() :
        four_momentum(mass, mass*rod_highgamma), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                rod_num_particles, rod_real_num,
                            comm_sptr, rod_length)
    {
        BOOST_TEST_MESSAGE("setup Round_rod bunch fixture highgamma");
        const double rod_radius = stdx * std::sqrt(2.0);
        bunch.set_sort_period(-1);
        MArray2d_ref local_particles(bunch.get_local_particles());
        // a ring of 8 particles around each longitudinal location
        int num_longitudinal = (rod_num_particles-1)/8;
        double dz = rod_length/(num_longitudinal-1);

        double r2o2 = std::sqrt(2.0)/2.0;

        double z = -rod_length/2.0;
        for (int i=4; i<rod_num_particles; i+=8, z+=dz) {
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

        // Add test particles at += stdx so mean and std are
        // not shifted
        local_particles[0][Bunch::x] = rod_probe;
        local_particles[0][Bunch::y] = 0.0;
        local_particles[0][Bunch::xp] = 0.0;
        local_particles[0][Bunch::yp] = 0.0;
        local_particles[0][Bunch::cdt] = 0.0;
        local_particles[0][Bunch::dpop] = 0.0;
        local_particles[0][Bunch::id] = 0.0;

        local_particles[1][Bunch::x] = -rod_probe;
        local_particles[1][Bunch::y] = 0.0;
        local_particles[1][Bunch::xp] = 0.0;
        local_particles[1][Bunch::yp] = 0.0;
        local_particles[1][Bunch::cdt] = 0.0;
        local_particles[1][Bunch::dpop] = 0.0;
        local_particles[1][Bunch::id] = 0.0;

        local_particles[2][Bunch::x] = 0.0;
        local_particles[2][Bunch::y] = rod_probe;
        local_particles[2][Bunch::xp] = 0.0;
        local_particles[2][Bunch::yp] = 0.0;
        local_particles[2][Bunch::cdt] = 0.0;
        local_particles[2][Bunch::dpop] = 0.0;
        local_particles[2][Bunch::id] = 0.0;

        local_particles[3][Bunch::x] = 0.0;
        local_particles[3][Bunch::y] = -rod_probe;
        local_particles[3][Bunch::xp] = 0.0;
        local_particles[3][Bunch::yp] = 0.0;
        local_particles[3][Bunch::cdt] = 0.0;
        local_particles[3][Bunch::dpop] = 0.0;
        local_particles[3][Bunch::id] = 0.0;
    }

    ~Round_rod_bunch_fixture_highgamma()
    {
        BOOST_TEST_MESSAGE("tear down Rod bunch fixture highgamma");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
};

BOOST_FIXTURE_TEST_CASE(apply_kv_round_highgamma, Round_rod_bunch_fixture_highgamma)
{
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    bunch.convert_to_state(Bunch::fixed_z_lab);
    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[2] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[2] << std::endl;
    double stdcdt = stds[4];
    logger << "Bunch stdcdt: " << stdcdt << std::endl;

    Space_charge_2d_kv space_charge;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    logger << "probe radius: " << probe_radius << std::endl;


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
    // volume factor is the charge contained in a cylinder of radius smaller than the probe
    double volume_factor = 1.0 - std::exp(-probe_radius*probe_radius/(2.0*stdx*stdx));
    logger << "volume_factor: " << volume_factor << std::endl;
    // fractional_charge includes the assumed longitudinal normal distribution
    double length_normalization = 1.0/L;
    logger << "length normalization factor: " << length_normalization << std::endl;
    double N = bunch.get_real_num()*volume_factor * length_normalization;

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(betagamma*betagamma*gamma)) *
            (step_length/bunch.get_local_particles()[0][Bunch::x]);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[0][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[1][Bunch::xp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[1][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[2][Bunch::yp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[2][Bunch::xp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[3][Bunch::yp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[3][Bunch::xp]), 1.0e-10);
}

BOOST_FIXTURE_TEST_CASE(apply_kv_round_gaussian_highgamma, Round_rod_bunch_fixture_highgamma)
{
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    bunch.convert_to_state(Bunch::fixed_z_lab);
    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[2] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[2] << std::endl;
    double stdcdt = stds[4];
    logger << "Bunch stdcdt: " << stdcdt << std::endl;

    Space_charge_2d_kv space_charge;
    // set to use gaussian charge density
    space_charge.set_longitudinal(Space_charge_2d_kv::longitudinal_gaussian);
    std::cout << "bassetti-erskine longitudinal: " << space_charge.get_longitudinal() << std::endl;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    logger << "probe radius: " << probe_radius << std::endl;


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
    // volume factor is the charge contained in a cylinder of radius smaller than the probe
    double volume_factor = 1.0 - std::exp(-probe_radius*probe_radius/(2.0*stdx*stdx));
    logger << "volume_factor: " << volume_factor << std::endl;
    // if fractional_charge uses longitudinal normal distribution
    double length_normalization = 1.0/(std::sqrt(2.0*mconstants::pi)*stdcdt*beta);

    logger << "length normalization factor: " << length_normalization << std::endl;
    double N = bunch.get_real_num()*volume_factor * length_normalization;

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(betagamma*betagamma*gamma)) *
            (step_length/bunch.get_local_particles()[0][Bunch::x]);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[0][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[1][Bunch::xp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[1][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[2][Bunch::yp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[2][Bunch::xp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[3][Bunch::yp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[3][Bunch::xp]), 1.0e-10);
}

BOOST_FIXTURE_TEST_CASE(apply_kv_round_offset_highgamma, Round_rod_bunch_fixture_highgamma)
{
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    bunch.convert_to_state(Bunch::fixed_z_lab);
    // add offset to bunch
    for (int i=0; i<bunch.get_local_num(); ++i) {
        bunch.get_local_particles()[i][Bunch::x] += 1.0;
        bunch.get_local_particles()[i][Bunch::y] += -2.0;
    }

    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[2] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[2] << std::endl;
    double stdcdt = stds[4];
    logger << "Bunch stdcdt: " << stdcdt << std::endl;

    Space_charge_2d_kv space_charge;
    std::cout << "bassetti-erskine longitudinal: " << space_charge.get_longitudinal() << std::endl;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    double probe_offset = bunch.get_local_particles()[0][Bunch::x]-mean[0];
    logger << "probe radius: " << probe_radius << std::endl;
    logger << "probe_offset: " << probe_offset << std::endl;

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
    // volume factor is the charge contained in a cylinder of radius smaller than the probe
    double volume_factor = 1.0 - std::exp(-probe_offset*probe_offset/(2.0*stdx*stdx));
    logger << "volume_factor: " << volume_factor << std::endl;
    // if fractional_charge includes the assumed longitudinal normal distribution
    // double length_normalization = 1.0/(std::sqrt(2.0*mconstants::pi)*stdcdt*beta);
    double length_normalization = 1.0/bunchlen;
    logger << "length normalization factor: " << length_normalization << std::endl;
    double N = bunch.get_real_num()*volume_factor * length_normalization;

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(betagamma*betagamma*gamma)) *
            (step_length/probe_offset);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[0][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[1][Bunch::xp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[1][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[2][Bunch::yp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[2][Bunch::xp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[3][Bunch::yp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[3][Bunch::xp]), 1.0e-10);
}

// I don't have a good test for elliptical distributions yet.  This is a place holder
struct Elliptical_rod_bunch_fixture_lowgamma
{
    Elliptical_rod_bunch_fixture_lowgamma() :
        four_momentum(mass, mass*rod_lowgamma), reference_particle(charge,
                four_momentum), comm_sptr(new Commxx), bunch(reference_particle,
                rod_num_particles, rod_real_num,
                            comm_sptr, rod_length)
    {
        BOOST_TEST_MESSAGE("setup Elliptical_rod bunch fixture lowgamma");
        const double rod_radius_x = stdx * std::sqrt(2.0)*1.01;
        const double rod_radius_y = stdx * std::sqrt(2.0)*0.99;
        bunch.set_sort_period(-1);
        MArray2d_ref local_particles(bunch.get_local_particles());
        // a ring of 8 particles around each longitudinal location
        int num_longitudinal = (rod_num_particles-1)/8;
        double dz = rod_length/(num_longitudinal-1);

        double r2o2 = std::sqrt(2.0)/2.0;

        double z = -rod_length/2.0;
        for (int i=4; i<rod_num_particles; i+=8, z+=dz) {
            local_particles[i][Bunch::x] = rod_radius_x;
            local_particles[i][Bunch::y] = 0.0;

            local_particles[i+1][Bunch::x] = rod_radius_x*r2o2;
            local_particles[i+1][Bunch::y] = rod_radius_y*r2o2;

            local_particles[i+2][Bunch::x] = 0.0;
            local_particles[i+2][Bunch::y] = rod_radius_y;

            local_particles[i+3][Bunch::x] = -rod_radius_x*r2o2;
            local_particles[i+3][Bunch::y] =  rod_radius_y*r2o2;

            local_particles[i+4][Bunch::x] = -rod_radius_x;
            local_particles[i+4][Bunch::y] = 0.0;

            local_particles[i+5][Bunch::x] = -rod_radius_x*r2o2;
            local_particles[i+5][Bunch::y] = -rod_radius_y*r2o2;

            local_particles[i+6][Bunch::x] = 0.0;
            local_particles[i+6][Bunch::y] = -rod_radius_y;

            local_particles[i+7][Bunch::x] = rod_radius_x*r2o2;
            local_particles[i+7][Bunch::y] = -rod_radius_y*r2o2;

            double rod_beta = reference_particle.get_beta();
            for (int j=i; j<i+8; ++j) {
                local_particles[j][Bunch::cdt] = z/rod_beta;
                local_particles[j][Bunch::xp] = 0.0;
                local_particles[j][Bunch::yp] = 0.0;
                local_particles[j][Bunch::dpop] = 0.0;
                local_particles[j][Bunch::id] = j;
            }
        }

        // Add test particles at += stdx so mean and std are
        // not shifted
        local_particles[0][Bunch::x] = stdx;
        local_particles[0][Bunch::y] = 0.0;
        local_particles[0][Bunch::xp] = 0.0;
        local_particles[0][Bunch::yp] = 0.0;
        local_particles[0][Bunch::cdt] = 0.0;
        local_particles[0][Bunch::dpop] = 0.0;
        local_particles[0][Bunch::id] = 0.0;

        local_particles[1][Bunch::x] = -stdx;
        local_particles[1][Bunch::y] = 0.0;
        local_particles[1][Bunch::xp] = 0.0;
        local_particles[1][Bunch::yp] = 0.0;
        local_particles[1][Bunch::cdt] = 0.0;
        local_particles[1][Bunch::dpop] = 0.0;
        local_particles[1][Bunch::id] = 0.0;

        local_particles[2][Bunch::x] = 0.0;
        local_particles[2][Bunch::y] = stdx;
        local_particles[2][Bunch::xp] = 0.0;
        local_particles[2][Bunch::yp] = 0.0;
        local_particles[2][Bunch::cdt] = 0.0;
        local_particles[2][Bunch::dpop] = 0.0;
        local_particles[2][Bunch::id] = 0.0;

        local_particles[3][Bunch::x] = 0.0;
        local_particles[3][Bunch::y] = -stdx;
        local_particles[3][Bunch::xp] = 0.0;
        local_particles[3][Bunch::yp] = 0.0;
        local_particles[3][Bunch::cdt] = 0.0;
        local_particles[3][Bunch::dpop] = 0.0;
        local_particles[3][Bunch::id] = 0.0;
    }

    ~Elliptical_rod_bunch_fixture_lowgamma()
    {
        BOOST_TEST_MESSAGE("tear down Rod bunch fixture lowgamma");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    unsigned long int seed;
};

// dummy test for when I come up with a good test for elliptical distributions
BOOST_FIXTURE_TEST_CASE(apply_kv_elliptical_offset_lowgamma, Elliptical_rod_bunch_fixture_lowgamma)
{
    return;
    const double time_fraction = 1.0;
    const double step_length = 0.1;
    const double beta = bunch.get_reference_particle().get_beta();
    const double betagamma = bunch.get_reference_particle().get_beta() * bunch.get_reference_particle().get_gamma();
    const double gamma = bunch.get_reference_particle().get_gamma();
    const double time_step = step_length/(beta*pconstants::c);
    const double bunchlen = bunch.get_z_period_length();

    Logger logger(0);

    bunch.convert_to_state(Bunch::fixed_z_lab);
    // add offset to bunch
    for (int i=0; i<bunch.get_local_num(); ++i) {
        bunch.get_local_particles()[i][Bunch::x] += 1.0;
        bunch.get_local_particles()[i][Bunch::y] += -2.0;
    }

    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d stds(Core_diagnostics::calculate_std(bunch, mean));
    logger << "Bunch means: " << mean[0] << ", " << mean[2] << std::endl;
    logger << "Bunch stds: " << stds[0] << ", " << stds[2] << std::endl;
    double stdcdt = stds[4];
    logger << "Bunch stdcdt: " << stdcdt << std::endl;

    Space_charge_2d_kv space_charge;
    std::cout << "bassetti-erskine longitudinal: " << space_charge.get_longitudinal() << std::endl;
    Step dummy_step(time_fraction);
    const int verbosity = 4;

    const double probe_radius = bunch.get_local_particles()[0][0];
    double probe_offset = bunch.get_local_particles()[0][Bunch::x]-mean[0];
    logger << "probe radius: " << probe_radius << std::endl;
    logger << "probe_offset: " << probe_offset << std::endl;

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
    // volume factor is the charge contained in a cylinder of radius smaller than the probe
    double volume_factor = 1.0 - std::exp(-probe_offset*probe_offset/(2.0*stdx*stdx));
    logger << "volume_factor: " << volume_factor << std::endl;
    // if fractional_charge includes the assumed longitudinal normal distribution
    // double length_normalization = 1.0/(std::sqrt(2.0*mconstants::pi)*stdcdt*beta);
    double length_normalization = 1.0/bunchlen;
    logger << "length normalization factor: " << length_normalization << std::endl;
    double N = bunch.get_real_num()*volume_factor * length_normalization;

    logger << "L: " << L << std::endl;
    logger << "N: " << N << std::endl;
    logger << "step_length: " << step_length << std::endl;
    logger << "betagamma: " << betagamma << std::endl;
    logger << "x: " << bunch.get_local_particles()[0][Bunch::x] << std::endl;
    double computed_dpop = ((2.0*N*pconstants::rp)/(betagamma*betagamma*gamma)) *
            (step_length/probe_offset);
    logger << "computed dpop: " << computed_dpop << std::endl;
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[0][Bunch::xp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[0][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[1][Bunch::xp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[1][Bunch::yp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[2][Bunch::yp], computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[2][Bunch::xp]), 1.0e-10);
    BOOST_CHECK_CLOSE(bunch.get_local_particles()[3][Bunch::yp], -computed_dpop, .01);
    BOOST_CHECK_LT(std::abs(bunch.get_local_particles()[3][Bunch::xp]), 1.0e-10);
}


// This checks that the bassetti-erskine calculation for nearly circular distributions is
// close to that of circular distributions.
BOOST_AUTO_TEST_CASE(close_to_circular)
{
    Space_charge_2d_kv sc1; // this one will be circular
    Space_charge_2d_kv sc2; // this one will be not quite circular

    const double stdx = 0.0005;
    const double stdcdt = 4.0;
    const double almost_circular = 1.0e-4;
    const double tolerance = 1.0e-4;
    const double really_close = 1.0e-13;

    bool isrnd1 = sc1.set_sigma(stdx, stdx, stdcdt);
    bool isrnd2 = sc2.set_sigma(stdx*(1.0+almost_circular), stdx*(1.0-almost_circular), stdcdt);

    BOOST_CHECK_EQUAL(isrnd1, 1);
    BOOST_CHECK_EQUAL(isrnd2, 0);

    // test near equality of electric field calculations over a grid
    const double minxy = -0.01;
    const double maxxy = 0.01;
    const int nsteps = 1001;
    const double dxy = (maxxy - minxy)/(nsteps-1);

    int cnt_tooclose = 0;

    double x = -minxy;
    for (int i=0; i<nsteps; ++i, x+=dxy) {
        double y = -minxy;
        for (int j=0; j<nsteps; ++j, y+=dxy) {
            double ex1, ey1;
            double ex2, ey2;
            sc1.normalized_efield(x, y, ex1, ey1);
            sc2.normalized_efield(x, y, ex2, ey2);

            // ex1 and ex2 should be close but not equal
            BOOST_CHECK(floating_point_equal(ex1, ex2, tolerance));
            BOOST_CHECK(floating_point_equal(ey1, ey2, tolerance));

            // count how many entries are really close
            if (floating_point_equal(ex1, ex2, really_close)) {
                ++cnt_tooclose;
            }
        }
    }
    BOOST_CHECK(cnt_tooclose < 100);
}

BOOST_AUTO_TEST_CASE(test_efield_calculation)
{
    // check normalized efield calculation against another program I wrote a long time ago
    // that implements the Bassetti-Erskine calculation
    const double stdcdt = 4.0;
    double x, y;
    double ex, ey;
    int i, j;
    const double sqrtpi = std::sqrt(pi);
    const double tolerance = 1.0e-6;

    Space_charge_2d_kv sc;

    double saved_ex[3][3] =
    {{-272.75132436615002, -282.47230821406555, -272.75132436615002},
     { 0.0000000000000000, 0.0, 0.0},
     { 272.75132436615002, 282.47230821406555, 272.75132436615002}};
    double saved_ey[3][3] =
    {{-667.98563127695525, 3.03641056817154449E-013, 667.98563127695525},
     {-677.47197461309258, 0.0,  677.47197461309258},
     {-667.98563127695525, 3.03641056817154449E-013, 667.98563127695525}};

    sc.set_sigma(stdx, stdy, stdcdt);
    for (i=0; i<3; ++i) {
        for (j=0; j<3; ++j) {
            x = -0.0001 + 0.0001*i;
            y = -0.0001 + 0.0001*j;
            sc.normalized_efield(x, y, ex, ey);
            BOOST_CHECK(floating_point_equal(ex*sqrtpi, saved_ex[i][j], tolerance));
            BOOST_CHECK(floating_point_equal(ey*sqrtpi, saved_ey[i][j], tolerance));
        }
    }
}
#endif
