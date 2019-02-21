#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice, MadX_reader
from synergia.bunch import Bunch
from synergia.utils import Commxx, Logger
from synergia.simulation import Bunch_simulator, Independent_stepper, Propagator, Lattice_simulator

import numpy as np

def read_nllens_lattice():
    ke = .0025
    energy = ke + pconstants.mp
    refpart = Reference_particle(1, pconstants.mp, energy)

    # read in the nllens lattice
    lattice = MadX_reader().get_lattice("channel", "./lattices/nllens_channel.seq")
    lattice.set_reference_particle(refpart)

    for elem in lattice.get_elements():
        elem.set_string_attribute("extractor_type", "libff")

    return lattice

def multipole_bunch(refpart):
    # create the bunch of particles matching the test particles
    # in the madx run
    macro_particles = 32
    real_particles = 1.0e10

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

    # repeat same pattern shifted 0.5mm down 0.25mm
    x2_offset = 0.0005
    y2_offset = -0.00025
    for i in range(16,32):
        lp[i, 0] = lp[i-16, 0] + x2_offset
        lp[i, 2] = lp[i-16, 2] + y2_offset

    return bunch

def run_nllens_channel():
    lattice = read_nllens_lattice()
    bunch = multipole_bunch(lattice.get_reference_particle())
    bunch_simulator = Bunch_simulator(bunch)
    stepper = Independent_stepper(lattice, 1, 1)
    propagator = Propagator(stepper)
    propagator.propagate(bunch_simulator, 1)
    m8p = np.load("./lattices/track_nllensone.npy")
    lp = bunch.get_local_particles()
    #numpart = lp.shape[0]
    numpart = bunch.get_local_num()
    assert(numpart == 32)
    for p in range(numpart):
        for j in range(4):
            assert_almost_equal(lp[p, j], m8p[p, j], 13)

def test_nllens():
    run_nllens_channel()

def test_nllens_map():
    lattice = read_nllens_lattice()
    # get knll and cnll parameters from first element
    nll = lattice.get_elements()[0]
    assert_equal(nll.get_type(), "nllens")
    knll = nll.get_double_attribute("knll")
    cnll = nll.get_double_attribute("cnll")
    k = 2*knll/cnll**2
    lattice_simulator = Lattice_simulator(lattice, 1)
    map = lattice_simulator.get_linear_one_turn_map()
    # the nllens should look like a quadrupole.  The
    # one turn map should be
    #  [ 1  0  0  0 ]
    #  [-k  1  0  0 ]
    #  [ 0  0  1  0 ]
    #  [ 0  0  k  0 ]
    # where k is 2*knll/cnll**2
    assert_almost_equal(-map[1,0], k, 15)
    assert_almost_equal(map[3,2], k, 15)
