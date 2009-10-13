#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "components/foundation/reference_particle.h"

const double tolerance = 1.0e-15;
const double mass = 100.0;
const double total_energy = 125.0;

BOOST_AUTO_TEST_CASE(construct)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    MArray1d state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        state[i] = 1.1 * i;
    }
    Reference_particle reference_particle(four_momentum, state);
}

BOOST_AUTO_TEST_CASE(get_four_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_four_momentum().get_total_energy(),
            four_momentum.get_total_energy(),tolerance);
}

BOOST_AUTO_TEST_CASE(get_state)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(total_energy);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],0.0,tolerance);
    }
}

BOOST_AUTO_TEST_CASE(get_beta)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_beta(),
            four_momentum.get_beta(),tolerance);
}

BOOST_AUTO_TEST_CASE(get_gamma)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_gamma(),
            four_momentum.get_gamma(),tolerance);
}

BOOST_AUTO_TEST_CASE(get_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_momentum(),
            four_momentum.get_momentum(),tolerance);
}

BOOST_AUTO_TEST_CASE(get_total_energy)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),
            four_momentum.get_total_energy(),tolerance);
}

BOOST_AUTO_TEST_CASE(set_four_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
    double new_total_energy = total_energy * 1.1;
    Four_momentum new_four_momentum(mass);
    new_four_momentum.set_total_energy(new_total_energy);
    reference_particle.set_four_momentum(new_four_momentum);
    BOOST_CHECK_CLOSE(reference_particle.get_four_momentum().get_total_energy(),
            new_total_energy,tolerance);
}

BOOST_AUTO_TEST_CASE(set_state)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    Reference_particle reference_particle(four_momentum);
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
    Reference_particle reference_particle(four_momentum);
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
    Reference_particle original_reference_particle(four_momentum, state);
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
    Reference_particle original_reference_particle(four_momentum);
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
