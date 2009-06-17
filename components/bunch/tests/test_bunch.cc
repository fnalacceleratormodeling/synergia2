#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "components/bunch/bunch.h"

const double tolerance = 1.0e-15;

const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 100;
const double real_num = 2.0e12;

BOOST_AUTO_TEST_CASE(construct)
{
	double units[] = {1.0,2.0,3.0,4.0,5.0,6.0};
    Reference_particle reference_particle(total_energy,units);
    Bunch bunch(reference_particle, proton_charge,
			total_num, real_num);
}
