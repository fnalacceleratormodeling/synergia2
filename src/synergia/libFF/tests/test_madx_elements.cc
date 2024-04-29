#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "catch2/matchers/catch_matchers.hpp"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

const double p0 = 1.5;
const double m = pconstants::mp;

const double default_tolerance = 1e-14;

struct propagator_fixture {
    Logger screen;
    Lattice lattice;
    std::unique_ptr<Propagator> propagator;
    std::unique_ptr<Bunch_simulator> sim;

    propagator_fixture(std::string const& file)
        : screen(0, LoggerV::INFO_TURN)
        , lattice(MadX_reader().get_lattice("machine",
                                            "./lattices/l_" + file + ".seq"))
        , propagator()
        , sim()
    {
        // Proton 1.5GeV/c
        auto fm_l = Four_momentum(pconstants::mp);
        fm_l.set_momentum(p0);
        auto ref_l = Reference_particle(1, fm_l);

        // Proton 1.5*1.05GeV/c
        auto fm_b = Four_momentum(pconstants::mp);
        fm_b.set_momentum(p0);
        auto ref_b = Reference_particle(1, fm_b);

        // set refpart to lattice
        lattice.set_reference_particle(ref_l);

        // cdt use drift (to agree with MadX)
        lattice.set_all_double_attribute("cdt_use_drift", 1.0);

        // propagator
        propagator = std::make_unique<Propagator>(
            lattice, Independent_stepper_elements(1));

        // bunch simulator
        sim = std::make_unique<Bunch_simulator>(
            Bunch_simulator::create_single_bunch_simulator(
                ref_b, 16, 1e10, Commxx()));
    }

    void
    propagate()
    {
        propagator->propagate(*sim, screen, 1);
    }

    Bunch&
    bunch()
    {
        return sim->get_bunch();
    }

    void
    print_lattice()
    {
        Logger l(0, LoggerV::DEBUG);
        lattice.print(l);
    }
};

karray2d_row
read_madx_output(std::string const& seq)
{
    karray2d_row out(seq, 16, 6);

    std::ifstream file("./lattices/madx_" + seq + ".out");
    std::string line;

    int row = 0;
    while (std::getline(file, line)) {
        if (row >= 26) {
            std::stringstream ss(line);

            double val;
            int col = 0;

            while (ss >> val) {
                if (col >= 2 && col <= 7) out(row - 26, col - 2) = val;
                ++col;
            }
        }
        ++row;
    }

    return out;
}

double
de_to_dp(double de)
{
    return sqrt(1 + de * de + 2 * de * sqrt(p0 * p0 + m * m) / p0) - 1.0;
}

void
propagate_libff(std::string const& seq, double tolerance = default_tolerance)
{
    std::cout << "libff propagate " << seq << "\n";

    propagator_fixture pf(seq);
    pf.print_lattice();

    // init particles
    auto& b = pf.bunch();

    double s2o2 = sqrt(2.0) / 2.0;
    double s3o2 = sqrt(3.0) / 2.0;
    double offset = 1.0e-3;

    double de = 0.05;
    double xp = 0.01;
    double yp = 0.015;
    double cdt = -0.01;

    // double dp = sqrt(1+de*de+2*de*sqrt(p0*p0+m*m)/p0) - 1.0;
    double dp = de_to_dp(de);

    b.checkout_particles();
    auto parts = b.get_host_particles();

    parts(0, 0) = offset;
    parts(0, 2) = 0.0;
    parts(0, 1) = xp;
    parts(0, 3) = yp;
    parts(0, 4) = cdt;
    parts(0, 5) = dp;

    parts(1, 0) = offset * s3o2;
    parts(1, 2) = offset * 0.5;
    parts(1, 1) = xp;
    parts(1, 3) = yp;
    parts(1, 4) = cdt;
    parts(1, 5) = dp;

    parts(2, 0) = offset * s2o2;
    parts(2, 2) = offset * s2o2;
    parts(2, 1) = xp;
    parts(2, 3) = yp;
    parts(2, 4) = cdt;
    parts(2, 5) = dp;

    parts(3, 0) = offset * 0.5;
    parts(3, 2) = offset * s3o2;
    parts(3, 1) = xp;
    parts(3, 3) = yp;
    parts(3, 4) = cdt;
    parts(3, 5) = dp;

    parts(4, 0) = 0.0;
    parts(4, 2) = offset;

    parts(5, 0) = -offset * 0.5;
    parts(5, 2) = offset * s3o2;

    parts(6, 0) = -offset * s2o2;
    parts(6, 2) = offset * s2o2;

    parts(7, 0) = -offset * s3o2;
    parts(7, 2) = offset * 0.5;

    parts(8, 0) = -offset;
    parts(8, 2) = 0.0;

    parts(9, 0) = -offset * s3o2;
    parts(9, 2) = -offset * 0.5;

    parts(10, 0) = -offset * s2o2;
    parts(10, 2) = -offset * s2o2;

    parts(11, 0) = -offset * 0.5;
    parts(11, 2) = -offset * s3o2;

    parts(12, 0) = 0.0;
    parts(12, 2) = -offset;

    parts(13, 0) = offset * 0.5;
    parts(13, 2) = -offset * s3o2;

    parts(14, 0) = offset * s2o2;
    parts(14, 2) = -offset * s2o2;

    parts(15, 0) = offset * s3o2;
    parts(15, 2) = -offset * 0.5;

    b.checkin_particles();

    // propagate
    pf.propagate();

    // print
    b.checkout_particles();
    parts = b.get_host_particles();

    // madx output
    auto madx = read_madx_output(seq);

    // check
    for (int p = 0; p < 16; ++p) {

        std::cout << std::setprecision(16);
        std::cout << std::scientific;

        std::cout << "libFF\tMadX\tmargin\n";

        for (int i = 0; i < 4; ++i) {
            std::cout << parts(p, i) << "\t" << madx(p, i) << "\t"
                      << fabs(parts(p, i) - madx(p, i)) << "\n";
        }

        std::cout << parts(p, 4) << "\t" << -madx(p, 4) << "\t"
                  << fabs(parts(p, 4) + madx(p, 4)) << "\n";

        std::cout << parts(p, 5) << "\t" << de_to_dp(madx(p, 5)) << "\t"
                  << fabs(parts(p, 5) - de_to_dp(madx(p, 5))) << "\n";

        std::cout << "\n";

        REQUIRE_THAT(parts(p, 0),
                     Catch::Matchers::WithinAbs(madx(p, 0), tolerance * 1e4));
        REQUIRE_THAT(parts(p, 1),
                     Catch::Matchers::WithinAbs(madx(p, 1), tolerance * 1e4));
        REQUIRE_THAT(parts(p, 2),
                     Catch::Matchers::WithinAbs(madx(p, 2), tolerance * 1e4));
        REQUIRE_THAT(parts(p, 3),
                     Catch::Matchers::WithinAbs(madx(p, 3), tolerance * 1e4));
        REQUIRE_THAT(parts(p, 4),
                     Catch::Matchers::WithinAbs(-madx(p, 4), tolerance * 1e6));
        REQUIRE_THAT(
            parts(p, 5),
            Catch::Matchers::WithinAbs(de_to_dp(madx(p, 5)), tolerance));
    }
}

TEST_CASE("drift", "[libFF][Elements]")
{
    propagate_libff("drift");
}

TEST_CASE("quad", "[libFF][Elements]")
{
    propagate_libff("quad");
}

TEST_CASE("sextupole", "[libFF][Elements]")
{
    propagate_libff("sext");
}

TEST_CASE("octupole", "[libFF][Elements]")
{
    propagate_libff("oct");
}

// long quad is tested with slightly reduced tolerance
// (1e-12 vs 1e-14)
TEST_CASE("quad_long", "[libFF][Elements]")
{
    propagate_libff("quad_long", 1e-12);
}

TEST_CASE("sextupole_long", "[libFF][Elements]")
{
    propagate_libff("sext_long");
}

TEST_CASE("octupole_long", "[libFF][Elements]")
{
    propagate_libff("oct_long");
}

TEST_CASE("kicker", "[libFF][Elements]")
{
    propagate_libff("kicker");
}

TEST_CASE("hkicker", "[libFF][Elements]")
{
    propagate_libff("hkicker");
}

TEST_CASE("vkicker", "[libFF][Elements]")
{
    propagate_libff("vkicker");
}

TEST_CASE("kicker_long", "[libFF][Elements]")
{
    propagate_libff("kicker_long");
}

TEST_CASE("hkicker_long", "[libFF][Elements]")
{
    propagate_libff("hkicker_long");
}

TEST_CASE("vkicker_long", "[libFF][Elements]")
{
    propagate_libff("vkicker_long");
}

TEST_CASE("rfc", "[libFF][Elements]")
{
    propagate_libff("rfc");
}

TEST_CASE("solenoid", "[libFF][Elements]")
{
    propagate_libff("solenoid");
}
