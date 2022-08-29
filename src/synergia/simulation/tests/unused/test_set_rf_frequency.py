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

def test_set_rf_frequency():
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "lattices/foborodobo128.madx")
    lattice_length = lattice.get_length()
    precision = -int(np.log(1.0e-13*lattice_length)/np.log(10.0))

    refpart = lattice.get_reference_particle()
    beta = refpart.get_beta()
    energy = refpart.get_total_energy()

    # propagate a particle.  With no kicks, 0 should propagate to 0
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)
    chef_beamline = lattice_simulator.get_chef_lattice().get_beamline()
    # clear the reference time for all elements
    total_reference_time = 0.0
    for ce in chef_beamline:
        total_reference_time = total_reference_time + ce.getReferenceTime()
        ce.setReferenceTime(0.0)
        if ce.Type() == "thinrfcavity":
            ce.setStrength(0.0)

    pr = beamline.Proton(energy)
    print "start: ", pr.State()
    chef_beamline.propagate(pr)
    print "end: ", pr.State()
    revolution_ctime = pr.get_cdt()
    print "revolution time from propagated particle: ", revolution_ctime
    print "lattice_simulator total reference ctime: ", total_reference_time

    inverse_t = synergia.foundation.pconstants.c/revolution_ctime
    print "frequency from propagation: ", inverse_t

    h = 1.0
    # naive_frequency is the frequency using just the length of
    # elements.
    naive_frequency = h * beta * synergia.foundation.pconstants.c/lattice_length
    print "naive frequency: ", naive_frequency

    assert_almost_equal(naive_frequency*1.0e-5, inverse_t*1.0e-5, precision)

    # creating the stepper will force registration setting
    # the RF frequency

    stepper_no_kicks = synergia.simulation.Independent_stepper(lattice, 1, 1)
    chef_beamline = stepper_no_kicks.get_lattice_simulator().get_chef_lattice().get_beamline()
    total_reference_time = 0.0
    for ce in chef_beamline:
        total_reference_time = total_reference_time + ce.getReferenceTime()
    print "stepper chef_beamline total reference time: ", total_reference_time

    # look in chef_beamline at the rf cavity
    frequency = None
    for chef_element in chef_beamline:
        if chef_element.Type() == "thinrfcavity":
            frequency = chef_element.getRadialFrequency()/(2.0 * np.pi)
            break
    if not frequency:
        raise RuntimeError, "Couldn't find thinrfcavity in chef_beamline"

    assert_almost_equal(frequency*1.0e-5, naive_frequency*1.0e-5, precision)

    # adding a kicker should change the RF frequency
    found_kicker = False
    for elem in lattice.get_elements():
        if elem.get_name() == "hc1":
            elem.set_double_attribute("kick", 0.02)
            found_kicker = True
            break

    if not found_kicker:
        raise RuntimeError("didn't find kicker")

    stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)
    # now the frequency should be different from the naive frequency
    chef_beamline = stepper.get_lattice_simulator().get_chef_lattice().get_beamline()

    # look in chef_beamline at the rf cavity
    frequency = None
    for chef_element in chef_beamline:
        if chef_element.Type() == "thinrfcavity":
            frequency = chef_element.getRadialFrequency()/(2.0 * np.pi)
            break

    assert_not_almost_equal(naive_frequency*1.0e-5, frequency*1.0e-5, precision)
