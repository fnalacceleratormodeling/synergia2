#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice, Lattice_element, chef_beamline_as_string, MadX_adaptor_map
from synergia.bunch import Bunch
from synergia.utils import Commxx, Logger
from synergia.simulation import Bunch_simulator, Independent_stepper, Propagator

import numpy as np

momentum = 1.5

macro_particles = 1
real_particles = 1.0e10

def get_the_bunch(refpart):
    commxx = Commxx()
    bunch = Bunch(refpart, macro_particles, real_particles, commxx)
    # particle with an angle in x.  With edges, there will be a y kick.
    # without edges, there will be no y kick
    lp = bunch.get_local_particles()
    lp[0,1] = 1.0e-3
    lp[0,2] = 1.0e-3
    return bunch

def run_bend(lattice):
    bunch = get_the_bunch(lattice.get_reference_particle())
    bunch_simulator = Bunch_simulator(bunch)
    stepper = Independent_stepper(lattice, 1, 1)
    propagator = Propagator(stepper)
    propagator.propagate(bunch_simulator, 1)
    return np.array(bunch.get_local_particles())

bend_length = 2.0
bend_angle = np.pi/12.0

# return reference particle of proton with momentum pc
def create_reference_particle(pc):
    four_momentum = Four_momentum(pconstants.mp)
    four_momentum.set_momentum(pc)
    reference_particle = Reference_particle(1, four_momentum)
    return reference_particle

def test_bend_with_edges():
    lattice = Lattice("foo", MadX_adaptor_map())
    bend = Lattice_element("sbend", "sb")
    bend.set_double_attribute("l", bend_length)
    bend.set_double_attribute("angle", bend_angle)
    bend.set_double_attribute("e1", bend_angle/2.0)
    bend.set_double_attribute("e2", bend_angle/2.0)
    zd = Lattice_element("drift", "zd")
    zd.set_double_attribute("l", 0.0)
    
    lattice.append(zd)
    lattice.append(bend)
    lattice.append(zd)
    lattice.set_reference_particle(create_reference_particle(momentum))
    particles = run_bend(lattice)
    print(particles[0,:])
    assert abs(particles[0, 3]) > 1.0e-10

def test_bend_no_edges():
    lattice = Lattice("foo", MadX_adaptor_map())
    bend = Lattice_element("sbend", "sb")
    bend.set_double_attribute("l", bend_length)
    bend.set_double_attribute("angle", bend_angle)
    bend.set_double_attribute("e1", bend_angle/2.0)
    bend.set_double_attribute("e2", bend_angle/2.0)
    bend.set_double_attribute("entry_edge_kick", 0.0)
    bend.set_double_attribute("exit_edge_kick", 0.0)
    zd = Lattice_element("drift", "zd")
    zd.set_double_attribute("l", 0.0)
    
    lattice.append(zd)
    lattice.append(bend)
    lattice.append(zd)
    
    lattice.set_reference_particle(create_reference_particle(momentum))
    particles = run_bend(lattice)
    print(particles[0,:])
    assert abs(particles[0,3]) < 1.0e-10
if __name__ == "__main__":
    print("test with edges")
    test_bend_with_edges()
    print()
    print("test no nedges")
    test_bend_no_edges()
