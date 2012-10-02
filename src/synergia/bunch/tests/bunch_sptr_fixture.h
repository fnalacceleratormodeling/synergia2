#ifndef BUNCH_SPTR_FIXTURE_H_
#define BUNCH_SPTR_FIXTURE_H_

#include "synergia/foundation/physical_constants.h"

const double mass = 100.0;
const double total_energy = 125.0;
const int local_num = 9;
int total_num;
const double real_num = 2.0e12;
const int turns = 17;
const double turn_length = 246.8;
const double partial_s = 123.4;

void
dummy_populate(Bunch &bunch)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; i += 1) {
            bunch.get_local_particles()[part][i] = 10.0 * part + (1.0 + part
                    * part / 1000.0) * i + 1000 * bunch.get_comm().get_rank();
        }
        bunch.get_local_particles()[part][Bunch::id] = part + total_num
                * bunch.get_comm().get_rank();
    }
}

struct Bunch_sptr_fixture
{
    Bunch_sptr_fixture() :
        comm_sptr(new Commxx),
        reference_particle(pconstants::electron_charge, mass,
                total_energy), bunch_sptr(new Bunch(reference_particle,
                local_num * comm_sptr->get_size(), real_num, comm_sptr))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        dummy_populate(*bunch_sptr);
        bunch_sptr->get_reference_particle().set_trajectory(turns, turn_length,
                partial_s);
    }
    ~Bunch_sptr_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Commxx_sptr comm_sptr;
    Reference_particle reference_particle;
    Bunch_sptr bunch_sptr;
};


#endif /* BUNCH_SPTR_FIXTURE_H_ */
