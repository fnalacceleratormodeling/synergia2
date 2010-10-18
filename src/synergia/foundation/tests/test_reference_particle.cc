#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "components/foundation/reference_particle.h"
#include "utils/xml_serialization.h"

const double tolerance = 1.0e-13;
const double mass = 100.0;
const double total_energy = 125.0;
const int charge = 1;
const double step_length = 1.234;
const int steps = 17;
const int turns = 7;

BOOST_AUTO_TEST_CASE(construct)
{
    Reference_particle reference_particle(charge, mass, total_energy);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
}

BOOST_AUTO_TEST_CASE(construct3)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    MArray1d state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        state[i] = 1.1 * i;
    }
    Reference_particle reference_particle(charge, four_momentum, state);
}

BOOST_AUTO_TEST_CASE(get_charge)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    BOOST_CHECK_EQUAL(reference_particle.get_charge(), charge);
}

BOOST_AUTO_TEST_CASE(get_four_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_four_momentum().get_total_energy(),
            four_momentum.get_total_energy(), tolerance);
}

BOOST_AUTO_TEST_CASE(get_state)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, total_energy);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i], 0.0, tolerance);
    }
}

BOOST_AUTO_TEST_CASE(get_beta)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_beta(),
            four_momentum.get_beta(), tolerance);
}

BOOST_AUTO_TEST_CASE(get_gamma)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_gamma(),
            four_momentum.get_gamma(), tolerance);
}

BOOST_AUTO_TEST_CASE(get_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_momentum(),
            four_momentum.get_momentum(), tolerance);
}

BOOST_AUTO_TEST_CASE(get_total_energy)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),
            four_momentum.get_total_energy(), tolerance);
}

BOOST_AUTO_TEST_CASE(increment_trajectory)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    for (int step = 0; step < steps; ++step) {
        reference_particle.increment_trajectory(step_length);
    }
    BOOST_CHECK_CLOSE(reference_particle.get_s(), steps*step_length,
            tolerance);
    BOOST_CHECK_CLOSE(reference_particle.get_trajectory_length(), steps*step_length,
            tolerance);
}

BOOST_AUTO_TEST_CASE(start_repetition)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    reference_particle.start_repetition();
    BOOST_CHECK_EQUAL(reference_particle.get_repetition(), 0);
    BOOST_CHECK_EQUAL(reference_particle.get_s(), 0);
    for (int step = 0; step < steps; ++step) {
        reference_particle.increment_trajectory(step_length);
    }
    reference_particle.start_repetition();
    BOOST_CHECK_EQUAL(reference_particle.get_repetition(), 1);
    BOOST_CHECK_EQUAL(reference_particle.get_s(), 0.0);
    BOOST_CHECK_CLOSE(reference_particle.get_repetition_length(),
            steps*step_length, tolerance);
}

BOOST_AUTO_TEST_CASE(set_trajectory)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    const double partial_s = 3 * step_length;
    reference_particle.set_trajectory(turns, steps * step_length, partial_s);
    BOOST_CHECK_EQUAL(reference_particle.get_repetition(), turns);
    BOOST_CHECK_CLOSE(reference_particle.get_repetition_length(),
            steps*step_length, tolerance);
    BOOST_CHECK_CLOSE(reference_particle.get_s(), partial_s, tolerance);
}

BOOST_AUTO_TEST_CASE(set_four_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    double new_total_energy = total_energy * 1.1;
    Four_momentum new_four_momentum(mass);
    new_four_momentum.set_total_energy(new_total_energy);
    reference_particle.set_four_momentum(new_four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_four_momentum().get_total_energy(),
            new_total_energy, tolerance);
}

BOOST_AUTO_TEST_CASE(set_state)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    MArray1d new_state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        new_state[i] = 1.1 * i;
    }
    reference_particle.set_state(new_state);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],new_state[i],
                tolerance);
    }
}

BOOST_AUTO_TEST_CASE(set_total_energy)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    double new_total_energy = total_energy * 1.1;
    reference_particle.set_total_energy(new_total_energy);
    BOOST_CHECK_CLOSE(reference_particle.get_four_momentum().get_total_energy(),
            new_total_energy,tolerance);
}

BOOST_AUTO_TEST_CASE(copy)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    MArray1d state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        state[i] = 1.1 * i;
    }
    Reference_particle
            original_reference_particle(charge, four_momentum, state);
    Reference_particle reference_particle(original_reference_particle);
    BOOST_CHECK_CLOSE(reference_particle.get_four_momentum().get_total_energy(),total_energy,
            tolerance);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],state[i],
                tolerance);
    }
}

BOOST_AUTO_TEST_CASE(copy2)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle original_reference_particle(charge, four_momentum);
    Reference_particle reference_particle(original_reference_particle);
    double new_total_energy = total_energy * 1.1;
    Four_momentum new_four_momentum(reference_particle.get_four_momentum());
    new_four_momentum.set_total_energy(new_total_energy);
    original_reference_particle.set_four_momentum(new_four_momentum);
    MArray1d new_state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        new_state[i] = 1.1 * i;
    }
    original_reference_particle.set_state(new_state);
    BOOST_CHECK_CLOSE(original_reference_particle.get_four_momentum().get_total_energy(),
            new_total_energy,tolerance);
    BOOST_CHECK_CLOSE(reference_particle.get_four_momentum().get_total_energy(),
            total_energy,tolerance);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(original_reference_particle.get_state()[i],
                new_state[i],tolerance);
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],0.0,tolerance);
    }
}

BOOST_AUTO_TEST_CASE(get_trajectory_length)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    for (int turn = 0; turn < turns; ++turn) {
        reference_particle.start_repetition();
        for (int step = 0; step < steps; ++step) {
            reference_particle.increment_trajectory(step_length);
        }
    }
    BOOST_CHECK_CLOSE(reference_particle. get_trajectory_length(),
            turns * steps * step_length, tolerance);
}

const double equal_tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(equal)
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);
    MArray1d new_state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        new_state[i] = 1.1 * i;
    }
    reference_particle1.set_state(new_state);
    Reference_particle reference_particle2(charge, four_momentum);
    reference_particle2.set_state(new_state);
    BOOST_CHECK(reference_particle1.equal(reference_particle2, equal_tolerance));
}

BOOST_AUTO_TEST_CASE(equal_different_four_momentum)
{
    Four_momentum four_momentum1(mass, total_energy);
    Reference_particle reference_particle1(charge, four_momentum1);
    Four_momentum four_momentum2(mass, total_energy * 1.1);
    Reference_particle reference_particle2(charge, four_momentum2);
    BOOST_CHECK(!reference_particle1.equal(reference_particle2, equal_tolerance));
}

BOOST_AUTO_TEST_CASE(equal_different_state)
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);
    MArray1d state1(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        state1[i] = 1.1 * i;
    }
    reference_particle1.set_state(state1);
    Reference_particle reference_particle2(charge, four_momentum);
    MArray1d state2(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        state2[i] = state1[i];
    }
    state2[1] *= 1.1;
    reference_particle2.set_state(state2);
    BOOST_CHECK(!reference_particle1.equal(reference_particle2, equal_tolerance));
}

BOOST_AUTO_TEST_CASE(equal_different_charge)
{
    Four_momentum four_momentum(mass);
    Reference_particle reference_particle1(charge, four_momentum);
    MArray1d new_state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        new_state[i] = 1.1 * i;
    }
    int charge2 = -1;
    reference_particle1.set_state(new_state);
    Reference_particle reference_particle2(charge2, four_momentum);
    reference_particle2.set_state(new_state);
    BOOST_CHECK(!reference_particle1.equal(reference_particle2, equal_tolerance));
}

BOOST_AUTO_TEST_CASE(test_serialize)
{
    Reference_particle reference_particle(charge, mass, total_energy);
    xml_save<Reference_particle > (reference_particle, "reference_particle.xml");

    Reference_particle loaded;
    xml_load<Reference_particle > (loaded, "reference_particle.xml");
    BOOST_CHECK(reference_particle.equal(loaded, tolerance));
}
