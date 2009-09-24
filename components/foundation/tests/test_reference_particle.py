#!/usr/bin/env python

#class Reference_particle
#{
#private:
#    double total_energy;
#    boost::multi_array<double, 1 > state;
#public:
#    Reference_particle(double total_energy);
#    Reference_particle(double total_energy, boost::const_multi_array_ref<
#            double, 1 > state);
#
#    void
#    set_total_energy(double total_energy);
#    void
#    set_state(boost::const_multi_array_ref<double, 1 > state);
#
#    double
#    get_total_energy();
#    boost::multi_array_ref<double, 1 >
#    get_state();
#};
import sys
sys.path.append('..')

from pyfoundation import Reference_particle
import numpy
from nose.tools import *

energy = 3.0
energy2 = 5.1
state = numpy.array([1.1,2.2,3.3,4.4,5.5,6.6],numpy.float)

def test_construct():
    r = Reference_particle(energy)

#def test_construct2():
#    r = Reference_particle(energy,state)

def test_set_get_total_energy():
    r = Reference_particle(energy)
    r.set_total_energy(energy2)
    assert_almost_equal(r.get_total_energy(),energy2)
    
def test_set_get_state():
    r = Reference_particle(energy)
    r.set_state(state)
    new_state = r.get_state()
    for i in range(0,6):
        assert_almost_equal(new_state[i],state[i])
