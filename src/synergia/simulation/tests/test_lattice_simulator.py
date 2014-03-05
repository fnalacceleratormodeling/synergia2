#!/usr/bin/env python

import synergia
from synergia.simulation import Lattice_simulator
import math
from nose.tools import *
from beamline import *

name = "fobodobo"
charge = 1
mass = synergia.foundation.pconstants.mp
p = 7.9447872893040543119
total_energy = math.sqrt(p**2 + mass**2)
tolerance = 1.0e-12
map_order = 2
n_cells = 8

bend_angle = 2 * math.pi / (2 * n_cells);
focus = 7;
sepn = 10;
quad_length = 0.2;
quad_strength = 1.0 / (focus * quad_length);
pct = 0.4;
bend_length = pct * (sepn - quad_length);
drift_length = (sepn - quad_length - bend_length) / 2.0;

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
        d.set_double_attribute("k1", -quad_strength)
        b = synergia.lattice.Lattice_element("sbend", "b")
        b.set_double_attribute("l", bend_length)
        b.set_double_attribute("angle", bend_angle)

        for cell in range(0, n_cells):
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

def test_adjust_tunes():
    f = Fixture()
    map_order = 1
    lattice_simulator = Lattice_simulator(f.lattice, map_order)

    horizontal_correctors = []
    vertical_correctors = []
    for element in f.lattice.get_elements():
        if element.get_type() == 'quadrupole':
            if element.get_double_attribute('k1') > 0.0:
                horizontal_correctors.append(element)
            else:
                vertical_correctors.append(element)

    new_horizontal_tune = 0.69
    new_vertical_tune = 0.15
    tolerance = 1.0e-6
    lattice_simulator.adjust_tunes(new_horizontal_tune, new_vertical_tune,
            horizontal_correctors, vertical_correctors, tolerance)

class Kick_Fixture:
    n_cells=32
    bend_angle = synergia.foundation.mconstants.pi/n_cells
    def __init__(self):
        four_momentum = synergia.foundation.Four_momentum(mass, 2.0)
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
        d.set_double_attribute("k1", -quad_strength)
        b = synergia.lattice.Lattice_element("sbend", "b")
        b.set_double_attribute("l", bend_length)
        b.set_double_attribute("angle", bend_angle)
        k = synergia.lattice.Lattice_element("hkicker", "k")
        k.set_double_attribute("kick", .0003)

        for cell in range(0, n_cells):
            self.lattice.append(f)
            self.lattice.append(o)
            self.lattice.append(b)
            self.lattice.append(o)
            self.lattice.append(d)
            self.lattice.append(o)
            self.lattice.append(b)
            self.lattice.append(o)
        self.lattice.append(k)

def test_get_closed_orbit():
    f = Kick_Fixture()
    map_order = 1
    lattice_simulator = Lattice_simulator(f.lattice, map_order)

    coords = lattice_simulator.get_closed_orbit()
    # with kick, the x-x' closed orbit better be nonzero.
    print "closed orbit:",
    #for i in range(6):
    #    print " ",coords[i],
    #print
    chef_beamline = lattice_simulator.get_chef_lattice().get_beamline()
    energy = chef_beamline.Energy()
    probe = Proton(energy)
    probe.set_x(coords[0])
    probe.set_npx(coords[1])
    probe.set_y(coords[2])
    probe.set_npy(coords[3])
    probe.set_cdt(coords[4])
    probe.set_ndp(coords[5])
    #print "probe before: ", probe.State()
    chef_beamline.propagate(probe)
    #print "probe after: ", probe.State()
    assert_almost_equal(coords[0], probe.get_x())
    assert_almost_equal(coords[1], probe.get_npx())




