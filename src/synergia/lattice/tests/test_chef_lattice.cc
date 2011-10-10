#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/chef_lattice.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/foundation/math_constants.h"
#include <basic_toolkit/PhysicsConstants.h>

const std::string name("fodo");
const int charge = 1;
const double mass = PH_NORM_mp;
const double total_energy = 8.9;
const double quad_length = 0.2;
const double quad_strength = 0.07;
const double drift_length = 3.0;
const double rf_length = 2.0;
const double rf_freq = 37.7e6;
const double tolerance = 1.0e-12;

struct Fodo_fixture
{
    Fodo_fixture() :
        four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice_sptr->set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);

        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
    }
    ~Fodo_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
};

const double bend_length = 0.1;
const int n_cells = 10;

struct Fobodobo_sbend_fixture
{
    Fobodobo_sbend_fixture() :
        four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice_sptr->set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);

        double bend_angle = 2 * mconstants::pi / (2 * n_cells);
        Lattice_element b("sbend", "b");
        b.set_double_attribute("l", bend_length);
        b.set_double_attribute("angle", bend_angle);

        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
    }
    ~Fobodobo_sbend_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
};

struct Fobodobo_sbend_markers_fixture
{
    Fobodobo_sbend_markers_fixture() :
        four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice_sptr->set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);

        double bend_angle = 2 * mconstants::pi / (2 * n_cells);
        Lattice_element b("sbend", "b");
        b.set_double_attribute("l", bend_length);
        b.set_double_attribute("angle", bend_angle);

        Lattice_element m("marker", "marker");

        lattice_sptr->append(f);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(b);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(d);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(b);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
    }
    ~Fobodobo_sbend_markers_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
};

struct Fobodobo_rbend_fixture
{
    Fobodobo_rbend_fixture() :
        four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice_sptr->set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);

        double bend_angle = 2 * mconstants::pi / (2 * n_cells);
        Lattice_element b("rbend", "b");
        b.set_double_attribute("l", bend_length);
        b.set_double_attribute("angle", bend_angle);

        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        lattice_sptr->append(b);
        lattice_sptr->append(o);
    }
    ~Fobodobo_rbend_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
};

struct Fobodobo_rbend_markers_fixture
{
    Fobodobo_rbend_markers_fixture() :
        four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice_sptr->set_reference_particle(reference_particle);
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);

        double bend_angle = 2 * mconstants::pi / (2 * n_cells);
        Lattice_element b("rbend", "b");
        b.set_double_attribute("l", bend_length);
        b.set_double_attribute("angle", bend_angle);

        Lattice_element m("marker", "marker");

        lattice_sptr->append(f);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(b);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(d);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(b);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
    }
    ~Fobodobo_rbend_markers_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
};

struct Chef_lattice_tester
{
    Chef_lattice chef_lattice;
    Chef_lattice_tester(Lattice_sptr lattice_sptr) :
        chef_lattice(lattice_sptr)
    {
    }
    Chef_elements
    get_chef_elements_from_slice(Lattice_element_slice const& slice)
    {
        return chef_lattice.get_chef_elements_from_slice(slice);
    }
};

struct RF_cavity_fixture
{
    RF_cavity_fixture() :
        four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice_sptr->set_reference_particle(reference_particle);
        Lattice_element rf("rfcavity", "rf");
        rf.set_double_attribute("l", rf_length);
        rf.set_double_attribute("freq", rf_freq);
        lattice_sptr->append(rf);
    }
    ~RF_cavity_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
};

BOOST_FIXTURE_TEST_CASE(construct, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
}

void
check_zero_reference_particle(Reference_particle const& reference_particle,
        double tolerance)
{
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK(std::abs(reference_particle.get_state()[i]) < tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_beamline_sptr, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    BmlPtr beamline_sptr = chef_lattice.get_beamline_sptr();
    const double tolerance = 1.0e-14;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_beamline_sptr_sbends, Fobodobo_sbend_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    BmlPtr beamline_sptr = chef_lattice.get_beamline_sptr();
    propagate_reference_particle(lattice_sptr->get_reference_particle(),
            beamline_sptr);
    const double tolerance = 1.0e-13;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_beamline_sptr_sbends_markers,
        Fobodobo_sbend_markers_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    BmlPtr beamline_sptr = chef_lattice.get_beamline_sptr();
    const double tolerance = 1.0e-13;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_beamline_sptr_rbends, Fobodobo_rbend_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    BmlPtr beamline_sptr = chef_lattice.get_beamline_sptr();
    const double tolerance = 1.0e-13;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_beamline_sptr_rbends_markers,
        Fobodobo_rbend_markers_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    BmlPtr beamline_sptr = chef_lattice.get_beamline_sptr();
    const double tolerance = 1.0e-13;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(have_sliced_beamline, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    BOOST_CHECK(!chef_lattice.have_sliced_beamline());
}

BOOST_FIXTURE_TEST_CASE(get_sliced_beamline_sptr_no_construct, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    bool caught = false;
    try {
        BmlPtr beamline_sptr = chef_lattice.get_sliced_beamline_sptr();
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

Lattice_element_slices
slice_lattice(Lattice & lattice, int slices_per_element)
{
    Lattice_element_slices slices;
    for (Lattice_elements::iterator it = lattice.get_elements().begin(); it
            != lattice.get_elements().end(); ++it) {
        double length = (*it)->get_length();
        if (length == 0.0) {
            Lattice_element_slice_sptr slice(new Lattice_element_slice(*(*it)));
            slices.push_back(slice);
        } else {
            double step_length = length / slices_per_element;
            for (int i = 0; i < slices_per_element; ++i) {
                double left = i * step_length;
                double right = (i + 1) * step_length;
                Lattice_element_slice_sptr slice(
                        new Lattice_element_slice(*(*it), left, right));
                slices.push_back(slice);
            }
        }
    }
    return slices;
}

BOOST_FIXTURE_TEST_CASE(get_sliced_beamline_sptr, Fodo_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    const int slices_per_element = 3;
    chef_lattice.construct_sliced_beamline(
            slice_lattice(*lattice_sptr, slices_per_element));
    BmlPtr beamline_sptr = chef_lattice.get_sliced_beamline_sptr();

    const double tolerance = 1.0e-14;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_sliced_beamline_sptr_sbends, Fobodobo_sbend_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    const int slices_per_element = 3;
    chef_lattice.construct_sliced_beamline(
            slice_lattice(*lattice_sptr, slices_per_element));
    BmlPtr beamline_sptr = chef_lattice.get_sliced_beamline_sptr();
    propagate_reference_particle(lattice_sptr->get_reference_particle(),
            beamline_sptr);
    const double tolerance = 1.0e-13;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_sliced_beamline_sptr_sbends_markers,
        Fobodobo_sbend_markers_fixture)
{
    Chef_lattice chef_lattice(lattice_sptr);
    const int slices_per_element = 3;
    chef_lattice.construct_sliced_beamline(
            slice_lattice(*lattice_sptr, slices_per_element));
    BmlPtr beamline_sptr = chef_lattice.get_sliced_beamline_sptr();
    const double tolerance = 1.0e-13;

    Reference_particle reference_particle(
            propagate_reference_particle(
                    lattice_sptr->get_reference_particle(), beamline_sptr));
    check_zero_reference_particle(reference_particle, tolerance);

    double synergia_length = lattice_sptr->get_length();
    double chef_length = beamline_sptr->OrbitLength(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()));
    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
}

//BOOST_FIXTURE_TEST_CASE(get_sliced_beamline_sptr_rbends, Fobodobo_rbend_fixture)
//{
//    Chef_lattice chef_lattice(lattice_sptr);
//    const int slices_per_element = 2;
//    chef_lattice.construct_sliced_beamline(
//            slice_lattice(*lattice_sptr, slices_per_element));
//    BmlPtr beamline_sptr = chef_lattice.get_sliced_beamline_sptr();
//    const double tolerance = 1.0e-13;
//
//    Reference_particle reference_particle(
//            propagate_reference_particle(
//                    lattice_sptr->get_reference_particle(), beamline_sptr));
//    check_zero_reference_particle(reference_particle, tolerance);
//
//    double synergia_length = lattice_sptr->get_length();
//    double chef_length = beamline_sptr->OrbitLength(
//            reference_particle_to_chef_particle(
//                    lattice_sptr->get_reference_particle()));
//    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
//}

//BOOST_FIXTURE_TEST_CASE(get_sliced_beamline_sptr_rbends_markers,
//        Fobodobo_rbend_markers_fixture)
//{
//    Chef_lattice chef_lattice(lattice_sptr);
//    BmlPtr beamline_sptr = chef_lattice.get_sliced_beamline_sptr();
//    const double tolerance = 1.0e-13;
//
//    Reference_particle reference_particle(
//            propagate_reference_particle(
//                    lattice_sptr->get_reference_particle(), beamline_sptr));
//    check_zero_reference_particle(reference_particle, tolerance);
//
//    double synergia_length = lattice_sptr->get_length();
//    double chef_length = beamline_sptr->OrbitLength(
//            reference_particle_to_chef_particle(
//                    lattice_sptr->get_reference_particle()));
//    BOOST_CHECK_CLOSE(synergia_length, chef_length, tolerance);
//}

double
chef_elements_length(Chef_elements const & chef_elements,
        Reference_particle const & reference_particle)
{
    double retval = 0.0;
    for (Chef_elements::const_iterator it = chef_elements.begin(); it
            != chef_elements.end(); ++it) {
        retval += (*it)->OrbitLength(
                reference_particle_to_chef_particle(reference_particle));
    }
    return retval;
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice1, Fodo_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element quad(*(lattice_sptr->get_elements().front()));
    Lattice_element_slice slice(quad);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(quad_length, chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice2, Fodo_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element quad(*(lattice_sptr->get_elements().front()));
    double begin = 0.0;
    double end = 1.0 / 3;
    Lattice_element_slice slice(quad, begin * quad_length, end * quad_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(quad_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice3, Fodo_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element quad(*(lattice_sptr->get_elements().front()));
    double begin = 1.0 / 3;
    double end = 2.0 / 3;
    Lattice_element_slice slice(quad, begin * quad_length, end * quad_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(quad_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice4, Fodo_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element quad(*(lattice_sptr->get_elements().front()));
    double begin = 2.0 / 3;
    double end = 1.0;
    Lattice_element_slice slice(quad, begin * quad_length, end * quad_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(quad_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice_compound1, RF_cavity_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element rfcavity(*(lattice_sptr->get_elements().front()));
    Lattice_element_slice slice(rfcavity);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(rf_length, chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
    BOOST_CHECK_EQUAL(chef_elements.size(), 3);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice_compound2, RF_cavity_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element rfcavity(*(lattice_sptr->get_elements().front()));
    double begin = 0.0;
    double end = 1.0 / 3;
    Lattice_element_slice slice(rfcavity, begin * rf_length, end * rf_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(rf_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
    BOOST_CHECK_EQUAL(chef_elements.size(), 1);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice_compound3, RF_cavity_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element rfcavity(*(lattice_sptr->get_elements().front()));
    double begin = 1.0 / 3;
    double end = 2.0 / 3;
    Lattice_element_slice slice(rfcavity, begin * rf_length, end * rf_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(rf_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
    BOOST_CHECK_EQUAL(chef_elements.size(), 3);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice_compound4, RF_cavity_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element rfcavity(*(lattice_sptr->get_elements().front()));
    double begin = 2.0 / 3;
    double end = 1.0;
    Lattice_element_slice slice(rfcavity, begin * rf_length, end * rf_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(rf_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
    BOOST_CHECK_EQUAL(chef_elements.size(), 1);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice_compound5, RF_cavity_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element rfcavity(*(lattice_sptr->get_elements().front()));
    double begin = 0.0;
    double end = 0.5;
    Lattice_element_slice slice(rfcavity, begin * rf_length, end * rf_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(rf_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
    BOOST_CHECK_EQUAL(chef_elements.size(), 2);
}

BOOST_FIXTURE_TEST_CASE(get_chef_elements_from_slice_compound6, RF_cavity_fixture)
{
    Chef_lattice_tester chef_lattice_tester(lattice_sptr);
    Lattice_element rfcavity(*(lattice_sptr->get_elements().front()));
    double begin = 0.5;
    double end = 1.0;
    Lattice_element_slice slice(rfcavity, begin * rf_length, end * rf_length);
    Chef_elements chef_elements(
            chef_lattice_tester.get_chef_elements_from_slice(slice));

    BOOST_CHECK_CLOSE(rf_length * (end - begin),
            chef_elements_length(chef_elements,
                    lattice_sptr->get_reference_particle()), tolerance);
    BOOST_CHECK_EQUAL(chef_elements.size(), 1);
}
