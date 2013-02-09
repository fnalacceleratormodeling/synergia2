#ifndef PROPAGATOR_FIXTURE_H_
#define PROPAGATOR_FIXTURE_H_

#include "lattice_fixture.h"

struct Propagator_fixture
{
    Propagator_fixture() :
            l(), space_charge(new Dummy_collective_operator("space_charge")), lattice_simulator(
                    l.lattice_sptr, 2), stepper_sptr(
                    new Split_operator_stepper(lattice_simulator, space_charge,
                            7)), propagator(stepper_sptr)
    {
    }
    ;

    ~Propagator_fixture()
    {
        BOOST_TEST_MESSAGE("teardown propagator fixture");
    }
    ;

    Lattice_fixture l;
    Dummy_collective_operator_sptr space_charge;
    Lattice_simulator lattice_simulator;
    Split_operator_stepper_sptr stepper_sptr;
    Propagator propagator;
};

#endif /* PROPAGATOR_FIXTURE_H_ */
