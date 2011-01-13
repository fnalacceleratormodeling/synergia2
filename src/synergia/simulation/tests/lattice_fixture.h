#ifndef LATTICE_FIXTURE_H_
#define LATTICE_FIXTURE_H_

#include "bunch_fixture.h"

const std::string name("foo");
const double quad_length = 0.2;
const double drift_length = 3.0;
const double bend_length = 4.0;

struct Lattice_fixture
{
    Lattice_fixture(): b(),lattice_sptr(new Lattice(name)) {
        BOOST_TEST_MESSAGE("setup lattice fixture");
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        Lattice_element m("monitor", "m");

        lattice_sptr->append(m);
        lattice_sptr->append(m);
        lattice_sptr->append(f);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(m);
        lattice_sptr->append(m);
        lattice_sptr->append(m);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(m);
        lattice_sptr->set_reference_particle(b.reference_particle);
    };

    ~Lattice_fixture(){
        BOOST_TEST_MESSAGE("teardown lattice fixture");

    };

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;
};

struct Lattice_fixture2
{
    Lattice_fixture2(): b(),lattice_sptr(new Lattice(name)) {
        BOOST_TEST_MESSAGE("setup lattice fixture");
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        Lattice_element m("monitor", "m");

        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        lattice_sptr->set_reference_particle(b.reference_particle);
    };

    ~Lattice_fixture2(){
        BOOST_TEST_MESSAGE("teardown lattice fixture");

    };

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;
};
#endif /* LATTICE_FIXTURE_H_ */
