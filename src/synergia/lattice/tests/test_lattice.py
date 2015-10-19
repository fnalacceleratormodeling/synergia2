#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from synergia.lattice import Lattice, Lattice_element
from synergia.foundation import Reference_particle
from nose.tools import *

name = "foo"
charge = -1
mass = 100.0
total_energy = 125.0
tolerance = 1.0e-12

def test_construct():
    lattice = Lattice(name)

def test_get_name():
    lattice = Lattice(name)
    assert_equal(lattice.get_name(), name)

def test_set_reference_particle():
    lattice = Lattice(name)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)

def test_has_reference_particle():
    lattice = Lattice(name)
    assert_equal(lattice.has_reference_particle(), False)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)
    assert_equal(lattice.has_reference_particle(), True)

def test_get_reference_particle():
    lattice = Lattice(name)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)
    assert_equal(lattice.get_reference_particle().equal(reference_particle, tolerance),
        True)
