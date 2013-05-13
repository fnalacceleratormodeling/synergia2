#!/usr/bin/env python

import sys
sys.path.append('..')
sys.path.append('../../convertors')

from foundation import Reference_particle, Four_momentum
import convertors

import numpy
from nose.tools import *

charge = 1
mass = 2.2;
total_energy = 3.0
total_energy2 = 5.1
four_momentum = Four_momentum(mass)
four_momentum.set_total_energy(total_energy)
state = numpy.array([1.1, 2.2, 3.3, 4.4, 5.5, 6.6], numpy.float)
step_length = 1.234;
steps = 17;
turns = 7;

def test_construct():
    r = Reference_particle(charge, mass, total_energy)

def test_construct2():
    r = Reference_particle(charge, four_momentum)

def test_construct3():
    r = Reference_particle(charge, four_momentum, state)

def test_get_charge():
    r = Reference_particle(charge, four_momentum)
    assert_equal(r.get_charge(), charge)

def test_set_get_four_momentum():
    r = Reference_particle(charge, four_momentum)
    r.get_four_momentum().set_total_energy(total_energy2)
    assert_almost_equal(r.get_four_momentum().get_total_energy(), total_energy2)

def test_set_get_four_momentum2():
    r = Reference_particle(charge, four_momentum)
    f2 = Four_momentum(mass)
    f2.set_total_energy(total_energy2)
    assert_almost_equal(r.get_four_momentum().get_total_energy(), total_energy)
    r.set_four_momentum(f2)
    assert_almost_equal(r.get_four_momentum().get_total_energy(), total_energy2)

def test_set_get_state():
    r = Reference_particle(charge, four_momentum)
    r.set_state(state)
    new_state = r.get_state()
    for i in range(0, 6):
        assert_almost_equal(new_state[i], state[i])

def test_set_get_total_energy():
    r = Reference_particle(charge, four_momentum)
    r.set_total_energy(total_energy2)
    assert_almost_equal(r.get_total_energy(), total_energy2)

def test_get_beta():
    r = Reference_particle(charge, four_momentum)
    assert_almost_equal(r.get_beta(), four_momentum.get_beta())

def test_get_gamma():
    r = Reference_particle(charge, four_momentum)
    assert_almost_equal(r.get_gamma(), four_momentum.get_gamma())

def test_get_momentum():
    r = Reference_particle(charge, four_momentum)
    assert_almost_equal(r.get_momentum(), four_momentum.get_momentum())

def test_set_get_trajectory():
    r = Reference_particle(charge, four_momentum)
    partial_s = 3 * step_length
    r.set_trajectory(turns, steps * step_length, partial_s)
    assert_equal(r.get_repetition(), turns)
    assert_almost_equal(r.get_repetition_length(), steps * step_length)
    assert_almost_equal(r.get_s_n(), partial_s)
    assert_almost_equal(r.get_s(),
                        turns * steps * step_length + partial_s)
