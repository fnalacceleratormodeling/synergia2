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

def read_solenoid_lattice():
    ke = .0025
    energy = ke + pconstants.mp
    refpart = Reference_particle(1, pconstants.mp, energy)

    # read in the nllens lattice
    lattice = MadX_reader().get_lattice("channel", "./lattices/solenoid_channel.seq")
    lattice.set_reference_particle(refpart)
    return lattice

def solenoid_bunch(refpart):
    # create the bunch of particles matching the test particles
    # in the madx run
    macro_particles = 16
    real_particles = 1.0e10

    commxx = Commxx()
    bunch = Bunch(refpart, macro_particles, real_particles, commxx)
    s2o2 = np.sqrt(2.0)/2.0
    s3o2 = np.sqrt(3.0)/2.0
    lp = bunch.get_local_particles()

    # 16 particles at radius of 1.0e-3 at 30degrees and 45 degrees offsets
    # around axis and give them transverse momentum so they spin
    offset = 1.0e-3
    lp[0, 0] = offset
    lp[0, 2] = 0.0
    lp[0, 1] = 0.0
    lp[0, 3] = offset/20.0
    
    lp[1, 0] = offset*s3o2
    lp[1, 2] = offset*0.5
    lp[1, 1] = -offset*0.5/20.0
    lp[1, 3] = offset*s3o2/20.0

    lp[2, 0] = offset*s2o2
    lp[2, 2] = offset*s2o2
    lp[2, 1] = -offset*s2o2/20.0
    lp[2, 3] = offset*s2o2/20.0

    lp[3, 0] = offset*0.5
    lp[3, 2] = offset*s3o2
    lp[3, 1] = -offset*s3o2/20.0
    lp[3, 3] = offset*0.5/20.0

    lp[4, 0] = 0.0
    lp[4, 2] = offset
    lp[4, 1] = -offset/20.0
    lp[4, 3] = 0.0

    lp[5, 0] = -offset*0.5
    lp[5, 2] =  offset*s3o2
    lp[5, 1] = -offset*s3o2/20.0
    lp[5, 3] = -offset*0.5/20.0

    lp[6, 0] = -offset*s2o2
    lp[6, 2] =  offset*s2o2
    lp[6, 1] = -offset*s2o2/20.0
    lp[6, 3] = -offset*s2o2/20.0

    lp[7, 0] = -offset*s3o2
    lp[7, 2] =  offset*0.5
    lp[7, 1] = -offset*0.5/20.0
    lp[7, 3] = -offset*s3o2/20.0

    lp[8, 0] = -offset
    lp[8, 2] = 0.0
    lp[8, 1] = 0.0
    lp[8, 3] = -offset/20.0

    lp[9, 0] = -offset*s3o2
    lp[9, 2] = -offset*0.5
    lp[9, 1] = offset*0.5/20.0
    lp[9, 3] = -offset*s3o2/20.0

    lp[10, 0] = -offset*s2o2
    lp[10, 2] = -offset*s2o2
    lp[10, 1] = offset*s2o2/20.0
    lp[10, 3] = -offset*s2o2/20.0

    lp[11, 0] = -offset*0.5
    lp[11, 2] = -offset*s3o2
    lp[11, 1] = offset*s3o2/20.0
    lp[11, 3] = -offset*0.5/20.0

    lp[12, 0] = 0.0
    lp[12, 2] = -offset
    lp[12, 1] = offset/20.0
    lp[12, 3] = 0.0

    lp[13, 0] =  offset*0.5
    lp[13, 2] = -offset*s3o2
    lp[13, 1] = offset*s3o2/20.0
    lp[13, 3] = offset*0.5/20.0

    lp[14, 0] =  offset*s2o2
    lp[14, 2] = -offset*s2o2
    lp[14, 1] = offset*s2o2/20.0
    lp[14, 3] = offset*s2o2/20.0

    lp[15, 0] =  offset*s3o2
    lp[15, 2] = -offset*0.5
    lp[15, 1] = offset*0.5/20.0
    lp[15, 3] = offset*s3o2/20.0

    return bunch

def run_solenoid_channel():
    print("foo bar")
    lattice = read_solenoid_lattice()
    bunch = solenoid_bunch(lattice.get_reference_particle())
    bunch_simulator = Bunch_simulator(bunch)
    stepper = Independent_stepper(lattice, 1, 1)
    print("before propagate: ", bunch.get_local_particles()[0,0:4])
    propagator = Propagator(stepper)
    propagator.propagate(bunch_simulator, 1)
    print("after propagate: ", bunch.get_local_particles()[0,0:4])
    m8p = np.load("./lattices/track_solenoidone.npy")
    lp = bunch.get_local_particles()
    numpart = lp.shape[0]
    assert(numpart == 16)
    for p in range(numpart):
        for j in range(4):
            assert_almost_equal(lp[p, j], m8p[p, j], 13)

def test_solenoid():
    run_solenoid_channel()

if __name__ == "__main__":
    run_solenoid_channel()
