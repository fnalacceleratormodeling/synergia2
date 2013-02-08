#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/utils/serialization.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "bunch_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_FIXTURE_TEST_CASE(construct, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
}

BOOST_FIXTURE_TEST_CASE(construct2, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Diagnostics_actions_sptr diagnostics_actions_sptr(new Diagnostics_actions);

    Bunch_simulator bunch_simulator(bunch_sptr, diagnostics_actions_sptr);
}

BOOST_FIXTURE_TEST_CASE(get_bunch, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    const double tolerance = 1.0e-12;
    BOOST_CHECK(
            bunch_simulator.get_bunch().get_reference_particle().get_four_momentum().equal( bunch_sptr->get_reference_particle().get_four_momentum(), tolerance));
}

BOOST_FIXTURE_TEST_CASE(get_bunch_sptr, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    const double tolerance = 1.0e-12;
    BOOST_CHECK(
            bunch_simulator.get_bunch_sptr()->get_reference_particle().get_four_momentum().equal( bunch_sptr->get_reference_particle().get_four_momentum(), tolerance));
}

BOOST_FIXTURE_TEST_CASE(get_diagnostics_actions, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    BOOST_CHECK_EQUAL(
            bunch_simulator.get_diagnostics_actions().get_bunch_sptr(),
            bunch_sptr);
}

BOOST_FIXTURE_TEST_CASE(get_diagnostics_actions_sptr, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    BOOST_CHECK_EQUAL(
            bunch_simulator.get_diagnostics_actions_sptr()->get_bunch_sptr(),
            bunch_sptr);
}

BOOST_FIXTURE_TEST_CASE(add_per_turn, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    bunch_simulator.add_per_turn(diagnostics_basic_sptr);
    const int period = 4;
    bunch_simulator.add_per_turn(diagnostics_basic_sptr, period);
}

BOOST_FIXTURE_TEST_CASE(add_per_turn2, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    bunch_simulator.add_per_turn(diagnostics_basic_sptr);
    std::list<int > prime_turns;
    prime_turns.push_back(2);
    prime_turns.push_back(3);
    prime_turns.push_back(5);

    bunch_simulator.add_per_turn(diagnostics_basic_sptr, prime_turns);
}

BOOST_FIXTURE_TEST_CASE(add_per_step, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    bunch_simulator.add_per_step(diagnostics_basic_sptr);
    const int period = 4;
    bunch_simulator.add_per_step(diagnostics_basic_sptr, period);
}

BOOST_FIXTURE_TEST_CASE(add_per_step2, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    bunch_simulator.add_per_step(diagnostics_basic_sptr);
    std::list<int > prime_steps;
    prime_steps.push_back(2);
    prime_steps.push_back(3);
    prime_steps.push_back(5);

    bunch_simulator.add_per_step(diagnostics_basic_sptr, prime_steps);
    const int period = 4;
    bunch_simulator.add_per_step(diagnostics_basic_sptr, prime_steps, period);
}

BOOST_FIXTURE_TEST_CASE(add_per_forced_diagnostics_step, Bunch_fixture)
{
    Bunch_sptr bunch_sptr(
            new Bunch(reference_particle, total_num, real_num, comm_sptr));
    Bunch_simulator bunch_simulator(bunch_sptr);
    const char filename[] = "dummy.h5";
    Diagnostics_basic_sptr diagnostics_basic_sptr(
            new Diagnostics_basic(filename));
    bunch_simulator.add_per_forced_diagnostics_step(diagnostics_basic_sptr);
    const int period = 4;
    bunch_simulator.add_per_forced_diagnostics_step(diagnostics_basic_sptr,
            period);
}

BOOST_FIXTURE_TEST_CASE(serialize, Bunch_fixture)
{
    {
        Bunch_sptr bunch_sptr(
                new Bunch(reference_particle, total_num, real_num, comm_sptr));
        Bunch_simulator bunch_simulator(bunch_sptr);
        xml_save<Bunch_simulator >(bunch_simulator, "bunch_simulator.xml");
    }

    {
        Bunch_simulator loaded;
        xml_load<Bunch_simulator >(loaded, "bunch_simulator.xml");
    }
}
