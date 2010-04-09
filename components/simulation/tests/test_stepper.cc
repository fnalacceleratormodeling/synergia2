#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/stepper.h"
#include "components/foundation/physical_constants.h"

const std::string name("foo");
const double mass = constants::mp;
const double total_energy = 125.0;
const double tolerance = 1.0e-12;
const double quad_length = 0.2;
const double drift_length = 3.0;
const double bend_length = 4.0;

BOOST_AUTO_TEST_CASE(construct)
{

    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(four_momentum);
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
    lattice.set_reference_particle(reference_particle);

    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice, 2);

    lattice.print();

    Split_operator_stepper stepper1(lattice_simulator, space_charge, 1);
    stepper1.print();

    Split_operator_stepper stepper2(lattice_simulator, space_charge, 2);
    stepper2.print();

    Split_operator_stepper stepper7(lattice_simulator, space_charge, 7);
    stepper7.print();

    Split_operator_stepper stepper10(lattice_simulator, space_charge, 10);

    stepper10.print();
}

