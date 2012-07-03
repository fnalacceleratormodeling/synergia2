#ifndef LATTICE_FIXTURE_H_
#define LATTICE_FIXTURE_H_

#include "bunch_fixture.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/serialization_files.h"
#include <cmath>

const std::string name("foo");
const double quad_length = 0.2;
const double quad_strength = 0.07;
const double drift_length = 3.0;
const double bend_length = 4.0;
const double drift_length_s = 0.1;

struct Lattice_fixture
{
    Lattice_fixture() :
        b(), lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup lattice fixture");
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        f.set_double_attribute("k1", quad_strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        d.set_double_attribute("k1", quad_strength);
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
    }
    ;

    ~Lattice_fixture()
    {
        BOOST_TEST_MESSAGE("teardown lattice fixture");

    }
    ;

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;
};

struct Lattice_fixture2
{
    Lattice_fixture2() :
        b(), lattice_sptr(new Lattice(name))
    {
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
    }
    ;

    ~Lattice_fixture2()
    {
        BOOST_TEST_MESSAGE("teardown lattice fixture");

    }
    ;

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;


};

struct Lattice_fixture3
{
    Lattice_fixture3() :
        b(), lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup lattice fixture");
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        Lattice_element m("monitor", "m");

        Lattice_element os("drift", "os");
        os.set_double_attribute("l", drift_length_s);


        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(os);
        lattice_sptr->append(os);
        lattice_sptr->append(d);
        lattice_sptr->append(os);
        lattice_sptr->append(o);
        lattice_sptr->set_reference_particle(b.reference_particle);
    }
    ;

    ~Lattice_fixture3()
    {
        BOOST_TEST_MESSAGE("teardown lattice fixture");

    }
    ;

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;


};
struct Lattice_fixture4
{
    Lattice_fixture4() :
        b(), lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup lattice fixture");
        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quad_length);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", drift_length);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quad_length);
        Lattice_element m("monitor", "m");

        Lattice_element os("drift", "os");
        os.set_double_attribute("l", drift_length_s);

        lattice_sptr->append(os);
        lattice_sptr->append(f);
        lattice_sptr->append(os);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        lattice_sptr->append(m);
        lattice_sptr->append(os);
        lattice_sptr->append(os);
        lattice_sptr->append(d);
        lattice_sptr->append(os);
        lattice_sptr->append(o);
        lattice_sptr->append(f);
        lattice_sptr->set_reference_particle(b.reference_particle);
    }
    ;

    ~Lattice_fixture4()
    {
        BOOST_TEST_MESSAGE("teardown lattice fixture");

    }
    ;

    Bunch_fixture b;
    Lattice_sptr lattice_sptr;


};



struct Fobodobo_sbend_fixture
{
    Fobodobo_sbend_fixture() :
                total_momentum(7.9447872893040543119),
                total_energy(
                        std::sqrt(
                                total_momentum * total_momentum
                                        + pconstants::mp * pconstants::mp)),
                four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name)), n_cells(8)
    {
        BOOST_TEST_MESSAGE("setup fixture");
        lattice_sptr->set_reference_particle(reference_particle);
        double bendangle = 2 * mconstants::pi / (2 * n_cells);
        double focus = 7;

        double sepn = 10;
        double quadlength = 0.2;
        double strength = 1.0 / (focus * quadlength);
        double pct = 0.4;
        double bendlength = pct * (sepn - quadlength);
        double driftlength = (sepn - quadlength - bendlength) / 2.0;

        Lattice_element f("quadrupole", "f");
        f.set_double_attribute("l", quadlength);
        f.set_double_attribute("k1", strength);
        Lattice_element o("drift", "o");
        o.set_double_attribute("l", driftlength);
        Lattice_element d("quadrupole", "d");
        d.set_double_attribute("l", quadlength);
        d.set_double_attribute("k1", -strength);

        Lattice_element b("sbend", "b");
        b.set_double_attribute("l", bendlength);
        b.set_double_attribute("angle", bendangle);

        for (int cell = 0; cell < n_cells; ++cell) {
            lattice_sptr->append(f);
            lattice_sptr->append(o);
            lattice_sptr->append(b);
            lattice_sptr->append(o);
            lattice_sptr->append(d);
            lattice_sptr->append(o);
            lattice_sptr->append(b);
            lattice_sptr->append(o);
        }
    }
    ~Fobodobo_sbend_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    double total_momentum;
    double total_energy;
    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
    int n_cells;
};

const int num_macro_particles = 4096;
const double num_real_particles = 1.0e11;
struct Foborodobo32_fixture
{
Foborodobo32_fixture() :
    lattice_sptr(new Lattice("foborodobo32")), comm_sptr(new Commxx),
    bunch_sptr()
  {
    BOOST_TEST_MESSAGE("setup Foborodobo_fixture");
    xml_load(*lattice_sptr, "foborodobo32_lattice.xml");
    bunch_sptr = Bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(),
				    num_macro_particles, num_real_particles, comm_sptr));
  }

  ~Foborodobo32_fixture()
  {
    BOOST_TEST_MESSAGE("teardown Foborodobo32 fixture");

  }

  Lattice_sptr lattice_sptr;
  Commxx_sptr comm_sptr;
  Bunch_sptr bunch_sptr;
};

#endif /* LATTICE_FIXTURE_H_ */
