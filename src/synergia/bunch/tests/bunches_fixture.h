#ifndef BUNCHES_FIXTURE_H_
#define BUNCHES_FIXTURE_H_

#include "synergia/bunch/tests/bunch_sptr_fixture.h"


struct Bunches_fixture
{
    Bunches_fixture() :
            num_bunches(3), reference_particle(pconstants::electron_charge,
                    mass, total_energy), bunches()
    {
        BOOST_TEST_MESSAGE("setup fixture");
	Commxx_sptr parent_sptr(new Commxx);
        Commxxs commxxs(generate_subcomms(parent_sptr, num_bunches));
        for (int i = 0; i < num_bunches; ++i) {
            Bunch_sptr bs(
                    new Bunch(reference_particle, local_num, real_num,
                            commxxs.at(i)));
            bunches.push_back(
                    bs);
            dummy_populate(*(bunches.back()));
        }
    }
  
    ~Bunches_fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    const int num_bunches;
    Reference_particle reference_particle;
    Bunches bunches;
};

#endif /* BUNCHES_FIXTURE_H_ */
