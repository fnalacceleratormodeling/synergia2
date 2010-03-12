#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/lattice/chef_lattice.h"
#include "components/lattice/chef_utils.h"
#include <basic_toolkit/PhysicsConstants.h>

const std::string name("fodo");
const double mass = PH_NORM_mp;
const double total_energy = 8.9;
const double quad_length = 0.2;
const double quad_strength = 0.07;
const double drift_length = 3.0;
const double tolerance = 1.0e-12;

struct Fodo_fixture
{
    Fodo_fixture() :
        four_momentum(mass, total_energy), reference_particle(four_momentum),
                lattice(name)
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice.set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);

        lattice.append(f);
        lattice.append(o);
        lattice.append(d);
        lattice.append(o);
    }
    ~Fodo_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Lattice lattice;
    Four_momentum four_momentum;
    Reference_particle reference_particle;
};

const double bend_length = 0.15;
const int n_cells = 3;
const double pi = 3.141592653589793;

struct Fobodobo_fixture
{
    Fobodobo_fixture() :
        four_momentum(mass, total_energy), reference_particle(four_momentum),
                lattice(name)
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice.set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);

        double bend_angle = 2 * pi / (2 * n_cells);
        Lattice_element b("sbend", "b");
        b.set_double_attribute("l", bend_length);
        b.set_double_attribute("angle", bend_angle);

        lattice.append(f);
        lattice.append(o);
        lattice.append(b);
        lattice.append(o);
        lattice.append(d);
        lattice.append(o);
        lattice.append(b);
        lattice.append(o);
    }
    ~Fobodobo_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Lattice lattice;
    Four_momentum four_momentum;
    Reference_particle reference_particle;
};

BOOST_FIXTURE_TEST_CASE(construct, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice);
}

BOOST_FIXTURE_TEST_CASE(construct2, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice,
            get_standard_lattice_element_to_chef_fn_map());
}

BOOST_FIXTURE_TEST_CASE(get_beamline_sptr, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice);
    BmlPtr beamline_sptr = chef_lattice.get_beamline_sptr();
    std::cout << "\nchef fodo\n";
    print_chef_beamline(beamline_sptr);
    // Not much of a test!!!
}

BOOST_FIXTURE_TEST_CASE(get_beamline_sptr_ring, Fobodobo_fixture)
{
    Chef_lattice chef_lattice(lattice);
    BmlPtr beamline_sptr = chef_lattice.get_beamline_sptr();
    std::cout << "\nchef fobodobo\n";
    print_chef_beamline(beamline_sptr);
    // Not much of a test!!!
}

