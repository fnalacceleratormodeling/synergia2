#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/chef_utils.h"

const int charge = 1;
const double mass = PH_NORM_mp;
const double total_energy = 8.9;

BOOST_AUTO_TEST_CASE(test_reference_particle_to_chef_jet_particle)
{
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    const int map_order = 3;
    JetParticle jet_particle(
            reference_particle_to_chef_jet_particle(reference_particle,
                    map_order));
    BOOST_CHECK_EQUAL(map_order, Jet__environment::getLastEnv()->maxWeight());
}

BOOST_AUTO_TEST_CASE(test_reference_particle_to_chef_jet_particle2)
{
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(charge, four_momentum);

    {
        const int map_order1 = 3;
        JetParticle jet_particle1(
                reference_particle_to_chef_jet_particle(reference_particle,
                        map_order1));
        BOOST_CHECK_EQUAL(map_order1,
                Jet__environment::getLastEnv()->maxWeight());
    }

    {
        const int map_order2 = 6;
        JetParticle jet_particle2(
                reference_particle_to_chef_jet_particle(reference_particle,
                        map_order2));
        BOOST_CHECK_EQUAL(map_order2,
                Jet__environment::getLastEnv()->maxWeight());
    }
}
