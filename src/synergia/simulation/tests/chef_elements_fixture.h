#ifndef CHEF_ELEMENTS_FIXTURE_H_
#define CHEF_ELEMENTS_FIXTURE_H_
#include <beamline/beamline_elements.h>
#include "synergia/simulation/chef_propagator.h"
#include "synergia/lattice/chef_utils.h"

const double drift_length = 1.2;
const double quad_length = 0.3;
const double quad_strength = 0.7;

struct Chef_elements_fixture
{
    Chef_elements_fixture()
    {
        BOOST_TEST_MESSAGE("setup Chef_elements fixture");
        ElmPtr drift1(new drift("drift1", drift_length));
        chef_elements.push_back(drift1);
        ElmPtr quad(new quadrupole("quad", quad_length, quad_strength));
        chef_elements.push_back(quad);
        ElmPtr drift2(new drift("drift2", drift_length));
        chef_elements.push_back(drift2);
    }
    ~Chef_elements_fixture()
    {
        BOOST_TEST_MESSAGE("teardown Chef_elements fixture");
    }
    Chef_elements chef_elements;
};

#endif /* CHEF_ELEMENTS_FIXTURE_H_ */
