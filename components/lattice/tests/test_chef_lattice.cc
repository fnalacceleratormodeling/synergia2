#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/lattice/chef_lattice.h"
#include <basic_toolkit/PhysicsConstants.h>

const std::string name("fodo");
const double mass = PH_NORM_mp;
const double total_energy = 8.9;
const double quad_length = 0.2;
const double drift_length = 3.0;
const double bend_length = 4.0;
const double tolerance = 1.0e-12;

struct Fixture
{
    Fixture() :
        four_momentum(mass, total_energy), reference_particle(four_momentum),
                lattice(name)
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice.set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);

        lattice.append(f);
        lattice.append(o);
        lattice.append(d);
        lattice.append(o);
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Lattice lattice;
    Four_momentum four_momentum;
    Reference_particle reference_particle;
};

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
    Chef_lattice chef_lattice(lattice);
}

BOOST_FIXTURE_TEST_CASE(construct2, Fixture)
{
    Chef_lattice chef_lattice(lattice,
            get_standard_lattice_element_to_chef_fn_map());
}

BOOST_FIXTURE_TEST_CASE(get_beamline_ptr, Fixture)
{
    Chef_lattice chef_lattice(lattice);
    BmlPtr beamline_ptr = chef_lattice.get_beamline_ptr();
    // Not much of a test!!!
}

