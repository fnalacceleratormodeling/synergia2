#!/usr/bin/env python

import synergia
import numpy as np
from synergia.foundation import pconstants
from nose.tools import *

# test that a tiny elliptical aperture cuts all the particles
def test_elliptical_aperture_cut_all():
    lattice = synergia.lattice.Lattice("testlattice")
    drift = synergia.lattice.Lattice_element("drift", "my_drift")
    drift.set_double_attribute("l", 1.0)
    drift.set_string_attribute("aperture_type", "elliptical")
    drift.set_double_attribute("elliptical_aperture_horizontal_radius", 0.001)
    drift.set_double_attribute("elliptical_aperture_vertical_radius", 0.0005)
    drift.set_string_attribute("extractor_type", "chef_propagate")
    lattice.append(drift)

    ref_momentum = 10.0
    pmom = synergia.foundation.Four_momentum(pconstants.mp)
    pmom.set_momentum(ref_momentum)
    refpart = synergia.foundation.Reference_particle(1, pmom)
    lattice.set_reference_particle(refpart)
    
    commxx = synergia.utils.Commxx()
    bunch = synergia.bunch.Bunch(refpart, 5, 1.0e10, commxx, 0.5)
    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
    
    local_particles = bunch.get_local_particles()
    local_num = bunch.get_local_num()

    for i in range(local_num):
        local_particles[i, 0] = 2.0

    stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)

    propagator = synergia.simulation.Propagator(stepper)
    propagator.propagate(bunch_simulator, 1, 1, 0)

    local_num = bunch.get_local_num()
    assert_equal(local_num, 0)

# test cutting the last particle
def test_elliptical_aperture_cut_last():
    lattice = synergia.lattice.Lattice("testlattice")
    drift = synergia.lattice.Lattice_element("drift", "my_drift")
    drift.set_double_attribute("l", 1.0)
    drift.set_string_attribute("aperture_type", "elliptical")
    drift.set_double_attribute("elliptical_aperture_horizontal_radius", 0.001)
    drift.set_double_attribute("elliptical_aperture_vertical_radius", 0.0005)
    drift.set_string_attribute("extractor_type", "chef_propagate")
    lattice.append(drift)

    ref_momentum = 10.0
    pmom = synergia.foundation.Four_momentum(pconstants.mp)
    pmom.set_momentum(ref_momentum)
    refpart = synergia.foundation.Reference_particle(1, pmom)
    lattice.set_reference_particle(refpart)
    
    commxx = synergia.utils.Commxx()
    bunch = synergia.bunch.Bunch(refpart, 5, 1.0e10, commxx, 0.5)
    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
    
    local_particles = bunch.get_local_particles()
    local_num = bunch.get_local_num()

    local_particles[local_num-1, 0] = 2.0

    stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)

    propagator = synergia.simulation.Propagator(stepper)
    propagator.propagate(bunch_simulator, 1, 1, 0)

    new_local_num = bunch.get_local_num()
    assert_equal(local_num-1, new_local_num)
    for i in range(new_local_num):
        print local_particles[i,:]

# test cutting all particles but last
def test_elliptical_aperture_cut_all_but_last():
    lattice = synergia.lattice.Lattice("testlattice")
    drift = synergia.lattice.Lattice_element("drift", "my_drift")
    drift.set_double_attribute("l", 1.0)
    drift.set_string_attribute("aperture_type", "elliptical")
    drift.set_double_attribute("elliptical_aperture_horizontal_radius", 0.001)
    drift.set_double_attribute("elliptical_aperture_vertical_radius", 0.0005)
    drift.set_string_attribute("extractor_type", "chef_propagate")
    lattice.append(drift)

    ref_momentum = 10.0
    pmom = synergia.foundation.Four_momentum(pconstants.mp)
    pmom.set_momentum(ref_momentum)
    refpart = synergia.foundation.Reference_particle(1, pmom)
    lattice.set_reference_particle(refpart)
    
    commxx = synergia.utils.Commxx()
    bunch = synergia.bunch.Bunch(refpart, 5, 1.0e10, commxx, 0.5)
    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
    
    local_particles = bunch.get_local_particles()
    local_num = bunch.get_local_num()

    for i in range(local_num-1):
        local_particles[i, 0] = 2.0

    orig_last_particle_id = int(local_particles[-1,6])

    stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)

    propagator = synergia.simulation.Propagator(stepper)
    propagator.propagate(bunch_simulator, 1, 1, 0)

    new_local_num = bunch.get_local_num()
    assert_equal(1, new_local_num)
    # remaining particle should be number the last of the original
    assert_equal(orig_last_particle_id, int(local_particles[0,6]))
