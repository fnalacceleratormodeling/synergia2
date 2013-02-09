#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/bunch_train_simulator.h"
#include "synergia/utils/serialization.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/tests/bunches_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double bunch_spacing = 1.7;

BOOST_FIXTURE_TEST_CASE(construct, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);
}

BOOST_FIXTURE_TEST_CASE(get_bunch_train, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);
    const double tolerance = 1.0e-12;
    for (int i = 0; i < bunch_train_sptr->get_size(); ++i) {
        BOOST_CHECK(
                bunch_train_simulator.get_bunch_train().get_bunches().at(i)->get_reference_particle().get_four_momentum().equal(
                        bunches.at(i)->get_reference_particle().get_four_momentum(), tolerance));
    }
}

BOOST_FIXTURE_TEST_CASE(get_bunch_train_sptr, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    BOOST_CHECK_EQUAL(bunch_train_simulator.get_bunch_train_sptr(), bunch_train_sptr);
}

BOOST_FIXTURE_TEST_CASE(get_diagnostics_actions, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    BOOST_CHECK_EQUAL(bunch_train_simulator.get_diagnostics_actionss().size(),
            bunches.size());
}

BOOST_FIXTURE_TEST_CASE(add_per_turn, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    const int which_bunch = 0;
    bunch_train_simulator.add_per_turn(which_bunch, diagnostics_basic_sptr);
    const int period = 4;
    bunch_train_simulator.add_per_turn(which_bunch, diagnostics_basic_sptr,
            period);
}

BOOST_FIXTURE_TEST_CASE(add_per_turn_bad_index, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    const int which_bunch = bunch_train_sptr->get_size();
    bool caught = false;
    try {
        bunch_train_simulator.add_per_turn(which_bunch, diagnostics_basic_sptr);
    }
    catch (std::out_of_range &) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_FIXTURE_TEST_CASE(add_per_turn2, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    const int which_bunch = 0;
    std::list<int > prime_turns;
    prime_turns.push_back(2);
    prime_turns.push_back(3);
    prime_turns.push_back(5);
    bunch_train_simulator.add_per_turn(which_bunch, diagnostics_basic_sptr,
            prime_turns);
}

BOOST_FIXTURE_TEST_CASE(add_per_step, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    const int which_bunch = 0;
    bunch_train_simulator.add_per_step(which_bunch, diagnostics_basic_sptr);
    const int period = 4;
    bunch_train_simulator.add_per_step(which_bunch, diagnostics_basic_sptr,
            period);
}

BOOST_FIXTURE_TEST_CASE(add_per_step_bad_index, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    const int which_bunch = bunch_train_sptr->get_size();
    bool caught = false;
    try {
        bunch_train_simulator.add_per_step(which_bunch, diagnostics_basic_sptr);
    }
    catch (std::out_of_range &) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_FIXTURE_TEST_CASE(add_per_step2, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    const int which_bunch = 0;
    std::list<int > prime_steps;
    prime_steps.push_back(2);
    prime_steps.push_back(3);
    prime_steps.push_back(5);
    bunch_train_simulator.add_per_step(which_bunch, diagnostics_basic_sptr,
            prime_steps);
}

BOOST_FIXTURE_TEST_CASE(add_per_forced_diagnostics_step, Bunches_fixture)
{
    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    const int which_bunch = 0;
    bunch_train_simulator.add_per_forced_diagnostics_step(which_bunch,
            diagnostics_basic_sptr);
    const int period = 4;
    bunch_train_simulator.add_per_forced_diagnostics_step(which_bunch,
            diagnostics_basic_sptr, period);
}

BOOST_FIXTURE_TEST_CASE(serialize, Bunches_fixture)
{
    {
        Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, bunch_spacing));
        Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

        xml_save<Bunch_train_simulator >(bunch_train_simulator, "bunch_train_simulator.xml");
    }

    {
        Bunch_train_simulator loaded;
        xml_load<Bunch_train_simulator >(loaded, "bunch_train_simulator.xml");
    }
}
