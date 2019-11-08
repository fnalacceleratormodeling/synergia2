#define CATCH_CONFIG_RUNNER
#include "synergia/utils/catch.hpp"

#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/independent_stepper_elements.h"


struct propagator_fixture
{
    Logger screen;
    Lattice lattice;
    Propagator propagator;
    Bunch_simulator sim;

    propagator_fixture(std::string const& seq)
        : screen(0, LoggerV::INFO_TURN)
        , lattice(MadX_reader().get_lattice(seq, "fodo.madx"))
        , propagator(lattice, Independent_stepper_elements(1))
        , sim( Bunch_simulator::create_single_bunch_simulator(
                    lattice.get_reference_particle(), 1, 1e09, Commxx()) )
    {
        auto ref = lattice.get_reference_particle();
        auto fm = ref.get_four_momentum();
        fm.set_momentum(fm.get_momentum()*0.95);
        //fm.set_momentum(3.0);
        ref.set_four_momentum(fm);

        //sim = Bunch_simulator::create_single_bunch_simulator(ref, 1, 1e09, Commxx());

        // propagate options (start turn, num turns)
        sim.set_turns(0, 1);
    }

    void propagate()
    { propagator.propagate(sim, screen); }

    Bunch& bunch()
    { return sim.get_bunch(); }
};



TEST_CASE("sbend", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_sbend");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}

TEST_CASE("sbend cf", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_cfsbend");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}

TEST_CASE("sbend cf 2", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_cfsbend2");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}

TEST_CASE("rfcavity", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_rfc");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}

#if 1
TEST_CASE("thin multipole 1 -- mad8 format", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_mp1");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}

TEST_CASE("thin multipole 2 -- mad8 format with tilts", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_mp2");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}
#endif


TEST_CASE("thin multipole 3 -- madx format", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_mp3");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}


TEST_CASE("thin multipole 4 -- madx format with tilts", "[libFF][Elements]")
{
    CHECK( true );

    propagator_fixture pf("seq_mp4");
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.1;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";
}


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    int result = Catch::Session().run(argc, argv);

    Kokkos::finalize();
    MPI_Finalize();

    return result;
}

