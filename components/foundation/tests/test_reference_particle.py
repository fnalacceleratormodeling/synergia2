#!/usr/bin/env python

import sys
sys.path.append('..')

from pyfoundation import Reference_particle, Four_momentum
import numpy
from nose.tools import *

mass = 2.2;
energy = 3.0
energy2 = 5.1
four_momentum = Four_momentum(mass)
four_momentum.set_total_energy(energy)
state = numpy.array([1.1,2.2,3.3,4.4,5.5,6.6],numpy.float)

def test_construct():
    r = Reference_particle(four_momentum)

def test_construct2():
    r = Reference_particle(four_momentum,state)

def test_set_get_four_momentum():
    r = Reference_particle(four_momentum)
    r.get_four_momentum().set_total_energy(energy2)
    assert_almost_equal(r.get_four_momentum().get_total_energy(),energy2)
    
def test_set_get_four_momentum2():
    r = Reference_particle(four_momentum)
    f2 = Four_momentum(mass)
    f2.set_total_energy(energy2)
    assert_almost_equal(r.get_four_momentum().get_total_energy(),energy)
    r.set_four_momentum(f2)
    assert_almost_equal(r.get_four_momentum().get_total_energy(),energy2)
    
def test_set_get_state():
    r = Reference_particle(four_momentum)
    r.set_state(state)
    new_state = r.get_state()
    for i in range(0,6):
        assert_almost_equal(new_state[i],state[i])
