#include "synergia/utils/catch.hpp"

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/utils/cereal_files.h"

const std::string name("foo");
const int charge = -1;
const double mass = 100.0;
const double total_energy = 125.0;
const double tolerance = 1.0e-12;

TEST_CASE("construct_lattice")
{
    Lattice lattice(name);
}

TEST_CASE("get_name_lattice")
{
    Lattice lattice(name);
    CHECK(lattice.get_name() == name);
}

TEST_CASE("set_reference_particle")
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
}

TEST_CASE("get_reference_particle")
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    CHECK(
        lattice.get_reference_particle().equal(reference_particle, tolerance));
}

const double quad_length = 0.2;
const double drift_length = 3.0;
const double bend_length = 4.0;

TEST_CASE("append_fodo")
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

    auto it = (lattice.get_elements()).begin();
    CHECK(it != lattice.get_elements().end());
    CHECK(it->get_name() == "f");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") ==
          Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") ==
          Approx(drift_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "d");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") ==
          Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") ==
          Approx(drift_length).margin(tolerance));
}

#if 0
TEST_CASE(derive_internal_attributes)
{
    Lattice_element b("rbend", "b");
    b.set_double_attribute("l", bend_length);
    double bend_angle = 2 * mconstants::pi / 24;
    b.set_double_attribute("angle", bend_angle);

    Lattice lattice(name);
    lattice.append(b);
    lattice.set_defaults();
    lattice.derive_internal_attributes();

    double arc_length = bend_angle * bend_length / (2
            * std::sin(bend_angle / 2));
    CHECK(lattice.get_length() == Approx(arc_length).margin(tolerance));
}
#endif

TEST_CASE("get_length")
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

    CHECK(lattice.get_length() ==
          Approx(2 * quad_length + 2 * drift_length).margin(tolerance));
}

TEST_CASE("get_total_angle1")
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

    CHECK(lattice.get_total_angle() == Approx(0.0).margin(tolerance));
}

TEST_CASE("get_total_angle2")
{
    const double pi = 3.141592653589793238462643;
    const int n_cells = 8;
    double bend_angle = 2 * pi / (2 * n_cells);

    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element b("sbend", "b");
    b.set_double_attribute("l", bend_length);
    b.set_double_attribute("angle", bend_angle);
    Lattice_element d("quadrupole", "d");
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

    CHECK(lattice.get_total_angle() == Approx(2 * pi).margin(tolerance));
}

TEST_CASE("lattice_from_lattice_element")
{
    Lattice_element f("drift", "foo");
    Lattice lattice(name);
    lattice.append(f);
    for (auto it = lattice.get_elements().begin();
         it != lattice.get_elements().end();
         ++it) {
        CHECK(it->has_lattice());
        CHECK(&(it->get_lattice()) == &lattice);
    }
}

TEST_CASE("copy_lattice")
{
    Lattice_element f("drift", "foo");
    const double foo_length = 1.0;
    f.set_double_attribute("l", foo_length);
    Lattice lattice(name);
    lattice.append(f);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);

    Lattice copied_lattice(lattice);

    CHECK(lattice.get_elements().size() ==
          copied_lattice.get_elements().size());
    CHECK(lattice.get_elements().begin()->get_name() ==
          copied_lattice.get_elements().begin()->get_name());

    CHECK(lattice.get_elements().begin()->get_length() ==
          Approx(foo_length).margin(tolerance));
    CHECK(copied_lattice.get_elements().begin()->get_length() ==
          Approx(foo_length).margin(tolerance));

    for (auto it = copied_lattice.get_elements().begin();
         it != copied_lattice.get_elements().end();
         ++it) {
        CHECK(it->has_lattice());
        CHECK(&(it->get_lattice()) == &copied_lattice);
    }

    const double new_length = 2.0;
    copied_lattice.get_elements().begin()->set_double_attribute("l",
                                                                new_length);

    CHECK(lattice.get_elements().begin()->get_length() ==
          Approx(foo_length).margin(tolerance));
    CHECK(copied_lattice.get_elements().begin()->get_length() ==
          Approx(new_length).margin(tolerance));

    const double new_energy = 2 * total_energy;
    Reference_particle new_reference_particle(charge, mass, new_energy);
    copied_lattice.set_reference_particle(new_reference_particle);
    CHECK(lattice.get_reference_particle().get_total_energy() ==
          Approx(total_energy).margin(tolerance));
    CHECK(copied_lattice.get_reference_particle().get_total_energy() ==
          Approx(new_energy).margin(tolerance));
}

TEST_CASE("copy_lattice2")
{
    Lattice lattice(name);
    Lattice copied_lattice(lattice);

    Reference_particle reference_particle(charge, mass, total_energy);
    copied_lattice.set_reference_particle(reference_particle);
}

TEST_CASE("copy_lattice_from_lattice_sptr")
{
    Lattice_element f("drift", "foo");
    const double foo_length = 1.0;
    f.set_double_attribute("l", foo_length);

    Lattice lattice(name);
    lattice.append(f);

    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);

    Lattice copied_lattice = lattice;

    CHECK(lattice.get_elements().size() ==
          copied_lattice.get_elements().size());
    CHECK(lattice.get_elements().begin()->get_name() ==
          copied_lattice.get_elements().begin()->get_name());

    CHECK(lattice.get_elements().begin()->get_length() ==
          Approx(foo_length).margin(tolerance));
    CHECK(copied_lattice.get_elements().begin()->get_length() ==
          Approx(foo_length).margin(tolerance));

    for (auto it = copied_lattice.get_elements().begin();
         it != copied_lattice.get_elements().end();
         ++it) {
        CHECK(it->has_lattice());
        CHECK(&(it->get_lattice()) == &copied_lattice);
    }
}

TEST_CASE("test_lsexpr")
{
#if 0
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
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);

    Lsexpr lsexpr(lattice.as_lsexpr());
    Lattice loaded(lsexpr);

    auto it = (loaded.get_elements()).begin();
    CHECK(it != loaded.get_elements().end());
    CHECK(it->get_name() == "f");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") == Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") == Approx(drift_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "d");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") == Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") == Approx(drift_length).margin(tolerance));
#endif
}

TEST_CASE("test_serialize1")
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
    json_save<Lattice>(lattice, "lattice1.json");

    Lattice loaded;
    json_load<Lattice>(loaded, "lattice1.json");

    auto it = loaded.get_elements().begin();
    CHECK(it != loaded.get_elements().end());
    CHECK(it->get_name() == "f");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") ==
          Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") ==
          Approx(drift_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "d");
    CHECK(it->get_type_name() == "quadrupole");
    CHECK(it->get_double_attribute("l") ==
          Approx(quad_length).margin(tolerance));
    ++it;

    CHECK(it->get_name() == "o");
    CHECK(it->get_type_name() == "drift");
    CHECK(it->get_double_attribute("l") ==
          Approx(drift_length).margin(tolerance));
}

TEST_CASE("test_serialize2")
{
    Lattice lattice(name);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice.set_reference_particle(reference_particle);
    json_save<Lattice>(lattice, "lattice2.json");

    Lattice loaded;
    json_load<Lattice>(loaded, "lattice2.json");

    CHECK(loaded.get_reference_particle().equal(reference_particle, tolerance));
}
