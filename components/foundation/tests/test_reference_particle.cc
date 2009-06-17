#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "reference_particle.h"

const double tolerance = 1.0e-15;

const double total_energy = 125.0;

BOOST_AUTO_TEST_CASE(construct)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle reference_particle(total_energy,units);
}

BOOST_AUTO_TEST_CASE(construct2)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
	double state[] = {1.1,2.2,3.2,4.3,5.5,6.6};
    Reference_particle reference_particle(total_energy,units,state);
}

BOOST_AUTO_TEST_CASE(get_total_energy)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle reference_particle(total_energy,units);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),total_energy,
        		tolerance);
}

BOOST_AUTO_TEST_CASE(get_state)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle reference_particle(total_energy,units);
    for(int i=0; i<6; ++i) {
    	BOOST_CHECK_CLOSE(reference_particle.get_state()[i],0.0,tolerance);
    }
}

BOOST_AUTO_TEST_CASE(get_units)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle reference_particle(total_energy,units);
    for(int i=0; i<6; ++i) {
    	BOOST_CHECK_CLOSE(reference_particle.get_units()[i],units[i],tolerance);
    }
}

BOOST_AUTO_TEST_CASE(set_total_energy)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle reference_particle(total_energy,units);
    double new_total_energy = total_energy * 1.1;
    reference_particle.set_total_energy(new_total_energy);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),new_total_energy,
        		tolerance);
}

BOOST_AUTO_TEST_CASE(set_state)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle reference_particle(total_energy,units);
    double new_state[] = {1.1,2.2,3.2,4.3,5.5,6.6};
    reference_particle.set_state(new_state);
    for(int i=0; i<6; ++i) {
    	BOOST_CHECK_CLOSE(reference_particle.get_state()[i],new_state[i],tolerance);
    }
}

BOOST_AUTO_TEST_CASE(copy)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle original_reference_particle(total_energy,units);
    Reference_particle reference_particle(original_reference_particle);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),total_energy,
        		tolerance);
    for(int i=0; i<6; ++i) {
    	BOOST_CHECK_CLOSE(reference_particle.get_state()[i],0.0,tolerance);
    	BOOST_CHECK_CLOSE(reference_particle.get_units()[i],units[i],tolerance);
    }
}

BOOST_AUTO_TEST_CASE(copy2)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle original_reference_particle(total_energy,units);
    Reference_particle reference_particle(original_reference_particle);
    double new_total_energy = total_energy * 1.1;
    original_reference_particle.set_total_energy(new_total_energy);
    double new_state[] = {1.1,2.2,3.2,4.3,5.5,6.6};
    original_reference_particle.set_state(new_state);
    BOOST_CHECK_CLOSE(original_reference_particle.get_total_energy(),
    		new_total_energy,tolerance);
    BOOST_CHECK_CLOSE(reference_particle.get_total_energy(),
    		total_energy,tolerance);
    for(int i=0; i<6; ++i) {
    	BOOST_CHECK_CLOSE(original_reference_particle.get_state()[i],
    			new_state[i],tolerance);
    	BOOST_CHECK_CLOSE(reference_particle.get_state()[i],0.0,tolerance);
    }
}
