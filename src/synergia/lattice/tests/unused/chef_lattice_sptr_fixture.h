#ifndef CHEF_LATTICE_FIXTURE_H_
#define CHEF_LATTICE_FIXTURE_H_

#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/chef_lattice_section.h"

const std::string name("fodo");
const double quad_length = 0.2;
const double quad_strength = 0.07;
const double drift_length = 3.0;

struct Fodo_fixture
{
    Fodo_fixture() :
        charge(1), mass(pconstants::mp), total_energy(8.9),
                four_momentum(mass, total_energy),
                reference_particle(charge, four_momentum),
                lattice_sptr(new Lattice(name))
    {
        BOOST_TEST_MESSAGE("setup Fodo_fixture");
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
        BOOST_TEST_MESSAGE("teardown Fodo_fixture");
    }

    int charge;
    double mass;
    double total_energy;
    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Lattice_sptr lattice_sptr;
};

struct Chef_lattice_sptr_fixture
{
    Chef_lattice_sptr_fixture() :
        fodo_fixture(),
                chef_lattice_sptr(new Chef_lattice(fodo_fixture.lattice_sptr))
    {
        BOOST_TEST_MESSAGE("setup Chef_lattice_fixture");
        for (Lattice_elements::iterator it =
                fodo_fixture.lattice_sptr->get_elements().begin(); it
                != fodo_fixture.lattice_sptr->get_elements().end(); ++it) {
            double length = (*it)->get_length();
            Lattice_element_slice_sptr slice1_sptr(
                    new Lattice_element_slice(*it, 0.0, 0.5 * length));
            Lattice_element_slice_sptr slice2_sptr(
                    new Lattice_element_slice(*it, 0.5 * length, length));
            slices.push_back(slice1_sptr);
            slices.push_back(slice2_sptr);
        }
        chef_lattice_sptr->construct_sliced_beamline(slices);
    }
    ~Chef_lattice_sptr_fixture()
    {
        BOOST_TEST_MESSAGE("teardown Chef_lattice_fixture");
    }

    Fodo_fixture fodo_fixture;
    Chef_lattice_sptr chef_lattice_sptr;
    Lattice_element_slices slices;
};

#endif /* CHEF_LATTICE_FIXTURE_H_ */
