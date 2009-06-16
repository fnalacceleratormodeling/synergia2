#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <iostream>

#include "four_momentum.h"

const double tolerance = 1.0e-15;

// some kinematics with simple decimal representations
const double mass = 100.0;
const double my_gamma = 1.25; // 5/4
const double beta = 0.6; // = sqrt(1-1/my_gamma^2) = 3/5
const double total_energy = 125.0;

BOOST_AUTO_TEST_CASE(construct)
{
    Four_momentum four_momentum(mass);
}

BOOST_AUTO_TEST_CASE(get_mass)
{
    Four_momentum four_momentum(mass);
    BOOST_CHECK_CLOSE(four_momentum.get_mass(),mass,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(get_total_energy)
{
    Four_momentum four_momentum(mass);
    BOOST_CHECK_CLOSE(four_momentum.get_total_energy(),mass,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(set_and_get_total_energy)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    BOOST_CHECK_CLOSE(four_momentum.get_total_energy(),total_energy,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(get_kinetic_energy)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    BOOST_CHECK_CLOSE(four_momentum.get_kinetic_energy(),
    		total_energy - mass,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(get_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    BOOST_CHECK_CLOSE(four_momentum.get_momentum(),
    		my_gamma * beta * mass,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(get_gamma)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    BOOST_CHECK_CLOSE(four_momentum.get_gamma(),
    		my_gamma,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(get_beta)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_total_energy(total_energy);
    BOOST_CHECK_CLOSE(four_momentum.get_beta(),
    		beta,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(set_kinetic_energy)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_kinetic_energy(total_energy-mass);
    BOOST_CHECK_CLOSE(four_momentum.get_total_energy(),total_energy,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(set_momentum)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_momentum(my_gamma*beta*mass);
    BOOST_CHECK_CLOSE(four_momentum.get_total_energy(),total_energy,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(set_gamma)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_gamma(my_gamma);
    BOOST_CHECK_CLOSE(four_momentum.get_total_energy(),total_energy,
    		tolerance);
}

BOOST_AUTO_TEST_CASE(set_beta)
{
    Four_momentum four_momentum(mass);
    four_momentum.set_beta(beta);
    BOOST_CHECK_CLOSE(four_momentum.get_total_energy(),total_energy,
    		tolerance);
}

//BOOST_AUTO_TEST_CASE(large_beta)
//{
//    Four_momentum four_momentum(mass);
//    double small = 1.0e-10;
//    double large_beta = 1.0 - small;
//    double large_gamma = 1.0/sqrt(2*small) + sqrt(2*small)/8.0; // Taylor expansion
//    four_momentum.set_beta(large_beta);
//    std::cout << "error is " << (four_momentum.get_gamma()-large_gamma)/large_gamma << std::endl;
//    BOOST_CHECK_CLOSE(four_momentum.get_gamma(),large_gamma,
//    		tolerance);
//}
//
