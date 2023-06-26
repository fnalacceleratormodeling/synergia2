#include "synergia/utils/catch.hpp"

#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/libFF/ff_element.h"


// Fixture that takes a string definition of an element,
// stuffs it into a lattice and creates a simulator/bunch ready to
// propagate, and a propagator to do the propagation.


struct propagator_fixture
{
    Logger screen;
    Lattice lattice;
    Propagator propagator;
    std::unique_ptr<Bunch_simulator> sim;

	
    propagator_fixture(std::string const& elem_def)
        : lattice(propagator_fixture::construct_test_lattice(elem_def))
        , screen(0, LoggerV::INFO_TURN)
        , propagator(Propagator(lattice, Independent_stepper_elements(1)))
        , sim()
    {
  
	    const Reference_particle& refpart = lattice.get_reference_particle();	    

        sim = std::make_unique<Bunch_simulator>(
                Bunch_simulator::create_single_bunch_simulator(
                    refpart, 8, 1e09, Commxx()));
    }

    Lattice
    construct_test_lattice(const std::string elem_def)
    {
        std::string lattice_head(R"foo(
beam, particle=proton, energy=0.8+pmass;
a: )foo");
        std::string lattice_tail(R"foo(
;
channel: sequence, refer=entry, l=0;
a, at=0.0;
endsequence;
)foo");

        // construct the lattice from the given element by sandwiching
        // it between the lattice_head and lattice_tail
	    MadX_reader reader;
	    reader.parse(lattice_head + elem_def + lattice_tail);
	    Lattice lattice(reader.get_lattice("channel"));

        return lattice;
    }

    // Propagate the same way it occurs on co_try
    void propagate()
    {
        for ( const auto& elm: propagator.get_lattice_elements() ) {
            FF_element::apply(elm, sim->get_bunch());
        }
    }

    Bunch& bunch()
    { return sim->get_bunch(); }

    void print_lattice()
    { Logger l(0, LoggerV::DEBUG); lattice.print(l); }
};

void propagate_test_elem(std::string const& elem_def, double tolerance)
{
    std::cout << "propagate test element " << elem_def << "\n";

    propagator_fixture pf(elem_def);
    pf.print_lattice();
    auto & b = pf.bunch();

    b.checkout_particles();
    auto parts = b.get_host_particles();
    for (int i=0; i<6; ++i) parts(0, i) = 0.0;
    b.checkin_particles();

    pf.propagate();

    b.checkout_particles();
    parts = b.get_host_particles();
    std::cout << std::scientific << std::setprecision(16);
    for (int i=0; i<6; ++i) std::cout << parts(0, i) << "\n";
    std::cout << "\n";

    // For propagating particle at 0, all the transverse elements should
    // remain at 0

    // on Ryzen 7,  CFsbends have a momentum of 4.5e-17
    for(int i=0; i<4; ++i) {
        CHECK (std::abs(parts(0, i)) < tolerance);
    }
}

TEST_CASE("sbend")
{
	propagate_test_elem("sbend, l=2, angle=pi/24", 2.0e-17);
}

TEST_CASE("CFsbend")
{
	propagate_test_elem("sbend, l=2, angle=pi/24, k1=1/16.2", 1.0e-11);
}

TEST_CASE("drift")
{
	propagate_test_elem("drift, l=10", 2.0e-17);
}

TEST_CASE("quadrupole")
{
	propagate_test_elem("quadrupole, l=4, k1=1/10", 2.0e-17);
}

TEST_CASE("sextupole")
{
	propagate_test_elem("sextupole, l=0.5, k2=0.25", 2.0e-17);
}

TEST_CASE("octupole")
{
	propagate_test_elem("octupole, l=0.5, k3=0.25", 2.0e-17);
}

TEST_CASE("skew-quadrupole")
{
	propagate_test_elem("quadrupole, l=4, k1=1/10, tilt=pi/4", 2.0e-17);
}

TEST_CASE("rfcavity")
{
	propagate_test_elem("rfcavity, l=2, volt=0.05", 2.0e-17);
}

