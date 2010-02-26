#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/lattice/lattice.h"

const std::string name("foo");
const double mass = 100.0;
const double total_energy = 125.0;

BOOST_AUTO_TEST_CASE(construct_lattice)
{
    Lattice lattice(name);
}

BOOST_AUTO_TEST_CASE(get_name_lattice)
{
    Lattice lattice(name);
    BOOST_CHECK_EQUAL(lattice.get_name(),name);
}

BOOST_AUTO_TEST_CASE(set_reference_particle)
{
    Lattice lattice(name);
    Reference_particle reference_particle(mass, total_energy);
    lattice.set_reference_particle(reference_particle);
}

BOOST_AUTO_TEST_CASE(has_reference_particle)
{
    Lattice lattice(name);
    BOOST_CHECK_EQUAL(lattice.has_reference_particle(),false);
    Reference_particle reference_particle(mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    BOOST_CHECK_EQUAL(lattice.has_reference_particle(),true);
}

BOOST_AUTO_TEST_CASE(get_reference_particle)
{
    Lattice lattice(name);
    Reference_particle reference_particle(mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    //    BOOST_CHECK_EQUAL(lattice.get_reference_particle(),reference_particle);
}
