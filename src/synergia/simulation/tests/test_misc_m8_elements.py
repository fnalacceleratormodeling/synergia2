#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice, Mad8_reader
from synergia.bunch import Bunch
from synergia.utils import Commxx, Logger
from synergia.simulation import Bunch_simulator, Independent_stepper, Propagator

import numpy as np

momentum = 1.5
multipole_sequence = "machine"
# read in a lattice, using the sequence name "machine"
# the reference particle for the lattice is a proton at momentum 1.5 GeV/c
def read_mad8_misc_element_lattice(elem_name):
    lattice = Mad8_reader().get_lattice(elem_name, "./lattices/misc_elements.lat")
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

def run_a_misc_element(elem_name):
    lattice = read_mad8_misc_element_lattice(elem_name)
    bunch = multipole_bunch(lattice.get_reference_particle())
    bunch_simulator = Bunch_simulator(bunch)
    stepper = Independent_stepper(lattice, 1, 1)
    propagator = Propagator(stepper)
    propagator.propagate(bunch_simulator, 1)
    m8p = np.load("./lattices/m8_"+elem_name+".npy")
    lp = bunch.get_local_particles()
    numpart = lp.shape[0]
    assert(numpart == 16)
    for p in range(numpart):
        for j in range(4):
            assert_almost_equal(lp[p, j], m8p[p, j], 10)
    
def test_base_quad():
    run_a_misc_element("m_base_quad")

def test_skew_quad():
    run_a_misc_element("m_skew_quad")

def test_tilt_quad():
    run_a_misc_element("m_tilt_quad")

