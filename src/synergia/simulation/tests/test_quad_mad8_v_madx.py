#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice, Lattice_element, Mad8_reader, MadX_reader
from synergia.bunch import Bunch
from synergia.utils import Commxx, Logger
from synergia.simulation import Bunch_simulator, Independent_stepper, Propagator

import numpy as np

momentum = 1.5
multipole_sequence = "machine"

# read in a lattice, using the sequence name "machine"
# the reference particle for the lattice is a proton at momentum 1.5 GeV/c
def get_m8_misc_lattice(elem_name):
    lattice = Mad8_reader().get_lattice(elem_name, "./lattices/misc_elements.lat")
    four_momentum = Four_momentum(pconstants.mp)
    four_momentum.set_momentum(momentum)
    refpart = Reference_particle(1, four_momentum)
    lattice.set_reference_particle(refpart)
    return lattice

def get_mx_misc_lattice(elem_name):
    lattice = MadX_reader().get_lattice(elem_name, "./lattices/misc_elements.seq")
    four_momentum = Four_momentum(pconstants.mp)
    four_momentum.set_momentum(momentum)
    refpart = Reference_particle(1, four_momentum)
    lattice.set_reference_particle(refpart)
    return lattice
    
def get_m8_multipole_lattice(elem_name):
    lattice = Mad8_reader().get_lattice("machine", "./lattices/" + elem_name + ".lat")
    four_momentum = Four_momentum(pconstants.mp)
    four_momentum.set_momentum(momentum)
    refpart = Reference_particle(1, four_momentum)
    lattice.set_reference_particle(refpart)
    return lattice

def get_mx_multipole_lattice(elem_name):
    lattice = MadX_reader().get_lattice("machine", "./lattices/" + elem_name + ".seq")
    four_momentum = Four_momentum(pconstants.mp)
    four_momentum.set_momentum(momentum)
    refpart = Reference_particle(1, four_momentum)
    lattice.set_reference_particle(refpart)
    return lattice


# create the bunch of particles matching the test particles
# in the madx run
macro_particles = 16
real_particles = 1.0e10

def multipole_bunch(refpart):
    commxx = Commxx()
    bunch = Bunch(refpart, macro_particles, real_particles, commxx)
    s2o2 = np.sqrt(2.0)/2.0
    s3o2 = np.sqrt(3.0)/2.0
    lp = bunch.get_local_particles()

    # 16 particles at radius of 1.0e-3 at 30degrees and 45 degrees offsets
    # around axis
    offset = 1.0e-3
    lp[0, 0] = offset
    lp[0, 2] = 0.0

    lp[1, 0] = offset*s3o2
    lp[1, 2] = offset*0.5

    lp[2, 0] = offset*s2o2
    lp[2, 2] = offset*s2o2

    lp[3, 0] = offset*0.5
    lp[3, 2] = offset*s3o2

    lp[4, 0] = 0.0
    lp[4, 2] = offset

    lp[5, 0] = -offset*0.5
    lp[5, 2] =  offset*s3o2

    lp[6, 0] = -offset*s2o2
    lp[6, 2] =  offset*s2o2

    lp[7, 0] = -offset*s3o2
    lp[7, 2] =  offset*0.5

    lp[8, 0] = -offset
    lp[8, 2] = 0.0

    lp[9, 0] = -offset*s3o2
    lp[9, 2] = -offset*0.5

    lp[10, 0] = -offset*s2o2
    lp[10, 2] = -offset*s2o2

    lp[11, 0] = -offset*0.5
    lp[11, 2] = -offset*s3o2

    lp[12, 0] = 0.0
    lp[12, 2] = -offset

    lp[13, 0] =  offset*0.5
    lp[13, 2] = -offset*s3o2

    lp[14, 0] =  offset*s2o2
    lp[14, 2] = -offset*s2o2

    lp[15, 0] =  offset*s3o2
    lp[15, 2] = -offset*0.5

    return bunch

def run_misc(elem_name):
    m8lattice = get_m8_misc_lattice(elem_name)
    mxlattice = get_mx_misc_lattice(elem_name)

    m8bunch = multipole_bunch(m8lattice.get_reference_particle())
    mxbunch = multipole_bunch(mxlattice.get_reference_particle())
    m8bunch_simulator = Bunch_simulator(m8bunch)
    mxbunch_simulator = Bunch_simulator(mxbunch)
    m8stepper = Independent_stepper(m8lattice, 1, 1)
    mxstepper = Independent_stepper(mxlattice, 1, 1)
    m8propagator = Propagator(m8stepper)
    mxpropagator = Propagator(mxstepper)
    m8propagator.propagate(m8bunch_simulator, 1)
    mxpropagator.propagate(mxbunch_simulator, 1)

    m8lp = m8bunch.get_local_particles()
    mxlp = mxbunch.get_local_particles()
    numpart = mxlp.shape[0]
    assert(numpart == 16)
    for p in range(numpart):
        for j in range(4):
            assert_almost_equal(m8lp[p, j], mxlp[p, j], 10)

def run_multipole(elem_name):
    m8lattice = get_m8_multipole_lattice(elem_name)
    mxlattice = get_mx_multipole_lattice(elem_name)

    m8bunch = multipole_bunch(m8lattice.get_reference_particle())
    mxbunch = multipole_bunch(mxlattice.get_reference_particle())
    m8bunch_simulator = Bunch_simulator(m8bunch)
    mxbunch_simulator = Bunch_simulator(mxbunch)
    m8stepper = Independent_stepper(m8lattice, 1, 1)
    mxstepper = Independent_stepper(mxlattice, 1, 1)
    m8propagator = Propagator(m8stepper)
    mxpropagator = Propagator(mxstepper)
    m8propagator.propagate(m8bunch_simulator, 1)
    mxpropagator.propagate(mxbunch_simulator, 1)

    m8lp = m8bunch.get_local_particles()
    mxlp = mxbunch.get_local_particles()
    numpart = mxlp.shape[0]
    assert(numpart == 16)
    for p in range(numpart):
        for j in range(4):
            assert_almost_equal(m8lp[p, j], mxlp[p, j], 10)
    
def test_base_quad():
    run_misc("m_base_quad")

def test_skew_quad():
    run_misc("m_skew_quad")

def test_tilt_quad():
    run_misc("m_tilt_quad")

def test_k1():
    run_multipole("mpole_k1")

def test_k1_tilt():
    run_multipole("mpole_k1_tilt")

def test_k2():
    run_multipole("mpole_k2")

def test_k2_tilt():
    run_multipole("mpole_k2_tilt")

def test_k3():
    run_multipole("mpole_k3")

def test_k3_tilt():
    run_multipole("mpole_k3_tilt")


def test_k4():
    run_multipole("mpole_k4")

def test_k4_tilt():
    run_multipole("mpole_k4_tilt")

def test_k5():
    run_multipole("mpole_k5")

def test_k5_tilt():
    run_multipole("mpole_k5_tilt")

def test_k6():
    run_multipole("mpole_k6")

def test_k6_tilt():
    run_multipole("mpole_k6_tilt")
