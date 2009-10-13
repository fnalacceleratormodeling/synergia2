#!/usr/bin/env python

import sys
sys.path.append('..')

from pyfoundation import Four_momentum
from nose.tools import *

mass = 3.0
total_energy = 27.4
value = 3.1415
beta_value = 0.678

def test_construct():
    f = Four_momentum(mass)

def test_construct2():
    f = Four_momentum(mass,total_energy)

def test_get_mass():
    f = Four_momentum(mass)
    assert_almost_equal(f.get_mass(),mass)

def test_set_get_total_energy():
    f = Four_momentum(mass)
    f.set_total_energy(value)
    assert_almost_equal(f.get_total_energy(),value)
    
def test_set_get_kinetic_energy():
    f = Four_momentum(mass)
    f.set_kinetic_energy(value)
    assert_almost_equal(f.get_kinetic_energy(),value)
    
def test_set_get_total_energy():
    f = Four_momentum(mass)
    f.set_total_energy(value)
    assert_almost_equal(f.get_total_energy(),value)
    
def test_set_get_beta():
    f = Four_momentum(mass)
    f.set_beta(beta_value)
    assert_almost_equal(f.get_beta(),beta_value)
    
