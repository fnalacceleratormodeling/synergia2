#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
import synergia
from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice, Mad8_reader, chef_beamline_as_string
from synergia.bunch import Bunch
from synergia.utils import Commxx, Logger
from synergia.simulation import Bunch_simulator, Independent_stepper, Propagator

import numpy as np
import beamline

def test_calculate_closed_orbit():
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "lattices/foborodobo128.madx")
    lattice_length = lattice.get_length()
    precision = -int(np.log(1.0e-13*lattice_length)/np.log(10.0))

    refpart = lattice.get_reference_particle()
    found_kicker = False
    for elem in lattice.get_elements():
        if elem.get_name() == "hc1":
            elem.set_double_attribute("kick", 0.02)
            found_kicker = True
            break

    if not found_kicker:
        raise RuntimeError("didn't find kicker")

    closed_orbit = synergia.simulation.calculate_closed_orbit(lattice)
    print("closed_orbit: ", closed_orbit)

    stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)

    propagator = synergia.simulation.Propagator(stepper)

    commxx = synergia.utils.Commxx()
    bunch = synergia.bunch.Bunch(lattice.get_reference_particle(), commxx.get_size(), 1.0e10, commxx)

    local_particles = bunch.get_local_particles()
    local_particles[0,0:6] = closed_orbit[:]

    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

    propagator.propagate(bunch_simulator, 1, 1, 0)

    for i in range(4):
        #print i, local_particles[0, i], closed_orbit[i]
        assert_almost_equal(local_particles[0, i], closed_orbit[i], precision)
