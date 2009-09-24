#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "reference_particle.h"

const double tolerance = 1.0e-15;

const double total_energy = 125.0;

BOOST_AUTO_TEST_CASE(construct)
{
    Reference_particle reference_particle(total_energy);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    boost::multi_array<double, 1 > state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        state[i] = 1.1 * i;
    }
    Reference_particle reference_particle(total_energy, state);
}

BOOST_AUTO_TEST_CASE(get_total_energy)
{
    Reference_particle reference_particle(total_energy);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),total_energy,
            tolerance);
}

BOOST_AUTO_TEST_CASE(get_state)
{
    Reference_particle reference_particle(total_energy);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],0.0,tolerance);
    }
}

BOOST_AUTO_TEST_CASE(set_total_energy)
{
    Reference_particle reference_particle(total_energy);
    double new_total_energy = total_energy * 1.1;
    reference_particle.set_total_energy(new_total_energy);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),new_total_energy,
            tolerance);
}

BOOST_AUTO_TEST_CASE(set_state)
{
    Reference_particle reference_particle(total_energy);
    boost::multi_array<double, 1 > new_state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        new_state[i] = 1.1 * i;
    }
    reference_particle.set_state(new_state);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],new_state[i],
                tolerance);
    }
}

BOOST_AUTO_TEST_CASE(copy)
{
    boost::multi_array<double, 1 > state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        state[i] = 1.1 * i;
    }
    Reference_particle original_reference_particle(total_energy, state);
    Reference_particle reference_particle(original_reference_particle);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),total_energy,
            tolerance);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],state[i],
                tolerance);
    }
}

BOOST_AUTO_TEST_CASE(copy2)
{
    Reference_particle original_reference_particle(total_energy);
    Reference_particle reference_particle(original_reference_particle);
    double new_total_energy = total_energy * 1.1;
    original_reference_particle.set_total_energy(new_total_energy);
    boost::multi_array<double, 1 > new_state(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        new_state[i] = 1.1 * i;
    }
    original_reference_particle.set_state(new_state);
    BOOST_CHECK_CLOSE(original_reference_particle.get_total_energy(),
            new_total_energy,tolerance);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),
            total_energy,tolerance);
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_CLOSE(original_reference_particle.get_state()[i],
                new_state[i],tolerance);
        BOOST_CHECK_CLOSE(reference_particle.get_state()[i],0.0,tolerance);
    }
}
