#!/usr/bin/env python

from synergia.lattice import Lattice, Lattice_element
from synergia.foundation import Reference_particle

name = "foo"
charge = -1
mass = 100.0
total_energy = 125.0
tolerance = 1.0e-12

def test_construct():
    lattice = Lattice(name)

def test_get_name():
    lattice = Lattice(name)
    assert(lattice.get_name() == name)

def test_set_reference_particle():
    lattice = Lattice(name)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)

def test_get_reference_particle():
    lattice = Lattice(name)
    reference_particle = Reference_particle(charge, mass, total_energy)
    lattice.set_reference_particle(reference_particle)
    assert(lattice.get_reference_particle().equal(reference_particle, tolerance))
