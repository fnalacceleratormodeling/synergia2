#include "synergia/utils/catch.hpp"

#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/independent_stepper_elements.h"


// Fixture that takes a string definition of an element,
// stuffs it into a lattice and creates a simulator/bunch ready to
// propagate, and a propagator to do the propagation.


struct propagator_fixture
{
    std::string test_elem;
    Logger screen;
    std::unique_ptr<Lattice> lattice;
    std::unique_ptr<Propagator> propagator;
    std::unique_ptr<Bunch_simulator> sim;

	static constexpr std::string lattice_head(R"foo(
beam, particle=proton, energy=0.8+pmass;
a: )foo");
	static constexpr std::string lattice_tail(R"foo(

channel: sequence, refer=entry, l=0;
a, at=0.0;
endsequence;
)foo");
	
    propagator_fixture(std::string const& elem_def)
        : test_elem(elem_def)
        , screen(0, LoggerV::INFO_TURN)
        , propagator()
        , lattice()
        , sim()
    {
        // construct the lattice from the given element
	    MadX_reader reader;
	    reader.parse(lattice_head + test_elem + lattice_tail);
	    lattice = reader.parse("channel");
  
	    Reference_particle refpart = lattice->get_reference_particle();
	    
	    propagator = Propagator(*lattice, Independent_stepper_elements(1));

        sim = std::make_unique<Bunch_simulator>(
                Bunch_simulator::create_single_bunch_simulator(
                    refpart, 8, 1e09, Commxx()));

    }

    void propagate()
    { *propagator.propagate(*sim, screen, 1); }

    Bunch& bunch()
    { return sim->get_bunch(); }

    void print_lattice()
    { Logger l(0, LoggerV::DEBUG); lattice->print(l); }
};

void propagate_libff(std::string const& seq)
{
    std::cout << "libff propagate " << seq << "\n";

    propagator_fixture pf(seq);
    pf.print_lattice();
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


TEST_CASE("sbend", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_sbend");
}

TEST_CASE("drift", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_drift");
}

TEST_CASE("quadrupole", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_quadrupole");
}

TEST_CASE("sbend cf", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_cfsbend");
}

TEST_CASE("sbend cf 2", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_cfsbend2");
}

TEST_CASE("octupole", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_octupole");
}

TEST_CASE("octupole with tilt", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_octupole2");
}

TEST_CASE("sextupole", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_sextupole");
}

TEST_CASE("sextupole with tilt", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_sextupole2");
}

TEST_CASE("rfcavity", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_rfc");
}

TEST_CASE("thin multipole 1 -- mad8 format", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_mp1");
}

TEST_CASE("thin multipole 2 -- mad8 format with tilts", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_mp2");
}

TEST_CASE("thin multipole 3 -- madx format", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_mp3");
}


TEST_CASE("thin multipole 4 -- madx format with tilts", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_mp4");
}

TEST_CASE("hkicker", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_hkicker");
}

TEST_CASE("vkicker", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_vkicker");
}

TEST_CASE("kicker", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_kicker");
}

TEST_CASE("long hkicker", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_long_hkicker");
}

TEST_CASE("long hkicker simple", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_long_hkicker_simple");
}

TEST_CASE("long vkicker", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_long_vkicker");
}

TEST_CASE("long vkicker simple", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_long_vkicker_simple");
}

TEST_CASE("long kicker", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_long_kicker");
}

TEST_CASE("long kicker simple", "[libFF][Elements]")
{
    CHECK( true );
    propagate_libff("seq_long_kicker_simple");
}



#if 0
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    int result = Catch::Session().run(argc, argv);

    Kokkos::finalize();
    MPI_Finalize();

    return result;
}

#endif
