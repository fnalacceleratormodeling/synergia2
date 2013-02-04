#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
from synergia.lattice import eliminate_markers, convert_monitors, convert_magnets, \
    combine_drifts, simplify_all
from synergia.lattice import Lattice, xml_load_lattice
from synergia.bunch import Bunch
from synergia.utils import Commxx, Logger
from synergia.simulation import Lattice_simulator, Independent_stepper

def fill(particles, delta):
    assert(particles.shape[0] == 7)
    for i in range(0, particles.shape[0]):
        particles[i, :] = 0.0
        particles[i, 6] = i
        if i > 0:
            particles[i, i - 1] = delta[i - 1]

def assert_particles_almost_equal(a, b, places):
    assert(a.shape[0] == b.shape[0])
    assert(a.shape[1] == b.shape[1])
    for i in range(0, a.shape[0]):
        for j in range(0, a.shape[1]):
            assert_almost_equal(a[i, j], b[i, j], places)

def general_exam(orig, simp, num_steps, map_order, places):
    total_num = 7
    real_num = 1.0
    bunch_orig = Bunch(orig.get_reference_particle(), total_num,
            real_num, Commxx());
    bunch_simp = Bunch(simp.get_reference_particle(), total_num,
            real_num, Commxx());
    delta = [0.001, 0.0001, 0.001, 0.0001, 0.1, 0.005]
    fill(bunch_orig.get_local_particles(), delta)
    fill(bunch_simp.get_local_particles(), delta)

    lattice_simulator_orig = Lattice_simulator(orig, map_order)
    lattice_simulator_simp = Lattice_simulator(simp, map_order)
    stepper_orig = Independent_stepper(lattice_simulator_orig, num_steps)
    stepper_simp = Independent_stepper(lattice_simulator_simp, num_steps)

    logger = Logger(0)
    verbosity = 0
    for step_orig, step_simp in zip(stepper_orig.get_steps(), stepper_simp.get_steps()):
        step_orig.apply(bunch_orig, verbosity, [], [], logger)
        step_simp.apply(bunch_simp, verbosity, [], [], logger)
        assert_particles_almost_equal(bunch_orig.get_local_particles(), \
                                      bunch_simp.get_local_particles(), places)

def test_booster():
    orig = Lattice()
    xml_load_lattice(orig, "lattices/fnal_booster.xml")
    simp = simplify_all(orig)

    general_exam(orig, simp, 96, 2, 13)


def test_debuncher():
    orig = Lattice()
    xml_load_lattice(orig, "lattices/fnal_debuncher.xml")
    simp = simplify_all(orig)

    general_exam(orig, simp, 204, 2, 13)

def test_main_injector():
    orig = Lattice()
    xml_load_lattice(orig, "lattices/fnal_main_injector.xml")
    simp = simplify_all(orig)

    general_exam(orig, simp, 416, 2, 11)

def test_main_injector_chef():
    orig = Lattice()
    xml_load_lattice(orig, "lattices/fnal_main_injector.xml")
    simp = simplify_all(orig)

    for element in orig.get_elements():
        element.set_string_attribute("extractor_type", "chef_propagate")
    for element in simp.get_elements():
        element.set_string_attribute("extractor_type", "chef_propagate")
    general_exam(orig, simp, 416, 1, 11)
