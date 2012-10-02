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
        Commxxs commxxs(generate_subcomms(num_bunches));
        for (int i = 0; i < num_bunches; ++i) {
            std::cout << "jfa1\n";
            std::cout << "jfa wtf: " << local_num << std::endl;
            Bunch_sptr bs(
                    new Bunch(reference_particle, local_num, real_num,
                            commxxs.at(i)));
            std::cout << "jfa1.1\n";
            bunches.push_back(
                    bs);
            std::cout << "jfa2\n";
            dummy_populate(*(bunches.back()));
        }
//        for (Commxxs::const_iterator it = commxxs.begin(); it != commxxs.end();
//                ++it) {
//            std::cout << "jfa1\n";
//            std::cout << "jfa wtf: " <<local_num << std::endl;
//            Bunch_sptr bs(new Bunch(reference_particle,
//                                    local_num * (*it)->get_size(), real_num,
//                                    Commxx_sptr(new Commxx)));
//            std::cout << "jfa1.1\n";
//            bunches.push_back(
//                    Bunch_sptr(
//                            new Bunch(reference_particle,
//                                    local_num * (*it)->get_size(), real_num,
//                                    *it)));
//            std::cout << "jfa2\n";
//            dummy_populate(*(bunches.back()));
//        }
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
