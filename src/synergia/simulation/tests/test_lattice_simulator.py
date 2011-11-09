#!/usr/bin/env python

import synergia
from synergia.simulation import Lattice_simulator
import math
from nose.tools import *

name = "fobodobo"
charge = 1
mass = synergia.foundation.pconstants.mp
total_energy = 8.9
quad_length = 0.2
quad_strength = 0.07
drift_length = 3.0
rf_length = 2.0
rf_freq = 37.7e6
tolerance = 1.0e-12
map_order = 2
n_cells = 8
quad_length = 0.2
quad_strength = 0.07
drift_length = 3.0
bend_length = 4.0

class Fixture:
    def __init__(self):
        four_momentum = synergia.foundation.Four_momentum(mass, total_energy)
        reference_particle = \
            synergia.foundation.Reference_particle(charge, four_momentum)
        self.lattice = synergia.lattice.Lattice(name)
        self.lattice.set_reference_particle(reference_particle)
        f = synergia.lattice.Lattice_element("quadrupole", "f")
        f.set_double_attribute("l", quad_length)
        f.set_double_attribute("k1", quad_strength)
        o = synergia.lattice.Lattice_element("drift", "o")
        o.set_double_attribute("l", drift_length)
        d = synergia.lattice.Lattice_element("quadrupole", "d")
        d.set_double_attribute("l", quad_length)
        d.set_double_attribute("k1", quad_strength)

        bend_angle = 2 * math.pi / (2 * n_cells)
        b = synergia.lattice.Lattice_element("sbend", "b")
        b.set_double_attribute("l", bend_length)
        b.set_double_attribute("angle", bend_angle)

        self.lattice.append(f)
        self.lattice.append(o)
        self.lattice.append(b)
        self.lattice.append(o)
        self.lattice.append(d)
        self.lattice.append(o)
        self.lattice.append(b)
        self.lattice.append(o)

def test_construct():
    f = Fixture()
    lattice_simulator = Lattice_simulator(f.lattice, map_order)