#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/lattice.h"

const std::string name("foo");
const int charge = -1;
const double mass = 100.0;
const double total_energy = 125.0;
const double tolerance = 1.0e-12;

#include "synergia/utils/xml_serialization.h"

BOOST_AUTO_TEST_CASE(construct_lattice)
{
    Lattice lattice(name);
}

BOOST_AUTO_TEST_CASE(get_name_lattice)
{
    Lattice lattice(name);
    BOOST_CHECK_EQUAL(lattice.get_name(), name);
}

BOOST_AUTO_TEST_CASE(set_reference_particle)
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
}

BOOST_AUTO_TEST_CASE(has_reference_particle)
{
    Lattice lattice(name);
    BOOST_CHECK_EQUAL(lattice.has_reference_particle(), false);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    BOOST_CHECK_EQUAL(lattice.has_reference_particle(), true);
}

BOOST_AUTO_TEST_CASE(get_reference_particle)
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    BOOST_CHECK(lattice.get_reference_particle().equal(reference_particle, tolerance));
}

const double quad_length = 0.2;
const double drift_length = 3.0;
const double bend_length = 4.0;

BOOST_AUTO_TEST_CASE(append_fodo)
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);

    Lattice_elements::const_iterator it = (lattice.get_elements()).begin();
    BOOST_CHECK(it != lattice.get_elements().end());
    BOOST_CHECK_EQUAL((*it)->get_name(), "f");
    BOOST_CHECK((*it)->get_type() == "quadrupole");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), quad_length, tolerance);
    ++it;

    BOOST_CHECK((*it)->get_name() == "o");
    BOOST_CHECK((*it)->get_type() == "drift");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), drift_length, tolerance);
    ++it;

    BOOST_CHECK((*it)->get_name() == "d");
    BOOST_CHECK((*it)->get_type() == "quadrupole");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), quad_length, tolerance);
    ++it;

    BOOST_CHECK((*it)->get_name() == "o");
    BOOST_CHECK((*it)->get_type() == "drift");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), drift_length, tolerance);
}

BOOST_AUTO_TEST_CASE(get_length)
{
    Lattice_element f("f", "quadrupole");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("o", "drift");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("d", "quadrupole");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);

    BOOST_CHECK_CLOSE(lattice.get_length(), 2*quad_length + 2*drift_length, tolerance);
}

BOOST_AUTO_TEST_CASE(get_total_angle1)
{
    Lattice_element f("f", "quadrupole");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("o", "drift");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("d", "quadrupole");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);

    BOOST_CHECK_CLOSE(lattice.get_total_angle(), 0.0, tolerance);
}

BOOST_AUTO_TEST_CASE(get_total_angle2)
{
    const double pi = 3.141592653589793238462643;
    const int n_cells = 8;
    double bend_angle = 2 * pi / (2 * n_cells);

    Lattice_element f("f", "quadrupole");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("o", "drift");
    o.set_double_attribute("l", drift_length);
    Lattice_element b("b", "sbend");
    b.set_double_attribute("l", bend_length);
    b.set_double_attribute("angle", bend_angle);
    Lattice_element d("d", "quadrupole");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    for (int i = 0; i < n_cells; ++i) {
        lattice.append(f);
        lattice.append(o);
        lattice.append(b);
        lattice.append(o);
        lattice.append(d);
        lattice.append(o);
        lattice.append(b);
        lattice.append(o);
    }

    BOOST_CHECK_CLOSE(lattice.get_total_angle(), 2*pi, tolerance);
}

BOOST_AUTO_TEST_CASE(test_serialize1)
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);
    xml_save<Lattice > (lattice, "lattice1.xml");

    Lattice loaded;
    xml_load<Lattice > (loaded, "lattice1.xml");

    BOOST_CHECK(!loaded.has_reference_particle());
    Lattice_elements::const_iterator it = (loaded.get_elements()).begin();
    BOOST_CHECK(it != loaded.get_elements().end());
    BOOST_CHECK_EQUAL((*it)->get_name(), "f");
    BOOST_CHECK((*it)->get_type() == "quadrupole");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), quad_length, tolerance);
    ++it;

    BOOST_CHECK((*it)->get_name() == "o");
    BOOST_CHECK((*it)->get_type() == "drift");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), drift_length, tolerance);
    ++it;

    BOOST_CHECK((*it)->get_name() == "d");
    BOOST_CHECK((*it)->get_type() == "quadrupole");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), quad_length, tolerance);
    ++it;

    BOOST_CHECK((*it)->get_name() == "o");
    BOOST_CHECK((*it)->get_type() == "drift");
    BOOST_CHECK_CLOSE((*it)->get_double_attribute("l"), drift_length, tolerance);
}

BOOST_AUTO_TEST_CASE(test_serialize2)
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    xml_save<Lattice > (lattice, "lattice2.xml");

    Lattice loaded;
    xml_load<Lattice > (loaded, "lattice2.xml");

    BOOST_CHECK(loaded.has_reference_particle());
    BOOST_CHECK(loaded.get_reference_particle().equal(reference_particle, tolerance));
}
