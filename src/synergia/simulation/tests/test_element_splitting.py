#!/usr/bin/env python

import sys

import pytest

from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice, MadX_reader
from synergia.bunch import Bunch
from synergia.utils import Commxx
from synergia.utils.parallel_utils import Logger, LoggerV
from synergia.simulation import Bunch_simulator, Independent_stepper_elements, Propagator

import numpy as np

mp = pconstants.mp
KE = 0.00250
places = 15

# Different lattices will have an element to test within a sequence named "machine"
# conditions from IOTA proton running
drift_txt = """
d: drift, l=2.0;
machine: sequence, l=2.0, refer=entry;
d, at = 0.0;
endsequence;
"""

m1r_txt = """
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
machine: sequence, l=0.3911403725, refer=entry;
m1r, at=0.0;
endsequence;
"""

m2r_txt = """
m2r: sbend,l:= 0.757985232,angle:= 1.047197551;
machine: sequence, l=0.757985232, refer=entry;
m2r, at=0.0;
endsequence;
"""

qa2r_txt = """
kq02 = 12.28401222;
kqa2r := kq02;
qa2r: quadrupole,l:= 0.21,k1:=kqa2r ;
machine: sequence, l=0.21, refer=entry;
qa2r, at=0.0;
endsequence;
"""

def read_iota_lattice_elem(e_txt):
    reader = MadX_reader()
    reader.parse(e_txt)
    lattice = reader.get_lattice('machine')

    refpart = Reference_particle(1, mp, mp+KE)
    lattice.set_reference_particle(refpart)
    return lattice


# create the bunch of particles matching the test particles
# in the madx run
macro_particles = 3*24 # 24 test particles in a ring at two different dp/p offsets
real_particles = 1.0e10

def make_sim(refpart):
    commxx = Commxx()
    sim = Bunch_simulator.create_single_bunch_simulator(
        refpart, macro_particles, real_particles, commxx)

    bunch = sim.get_bunch(0, 0)

    s2o2 = np.sqrt(2.0)/2.0
    s3o2 = np.sqrt(3.0)/2.0
    s15d = (np.sqrt(6) - np.sqrt(2))/4.0 # sin 15degrees
    c15d = (np.sqrt(6) + np.sqrt(2))/4.0 # cos 15degrees
    s75d = c15d
    c75d = s15d

    lp = bunch.get_particles_numpy()

    # 24 particles at radius of 1.0e-3 at 30degrees and 45 degrees offsets
    # around axis
    offset = 0.00213*3 # 3 times rms y
    dpopoffset = 2.1e-3 * 2.5 # 2.5 times rms dp/p spread

    lp[0, 0] = offset
    lp[0, 2] = 0.0

    lp[1, 0] = offset*c15d # 15 degrees
    lp[1, 2] = offset*s15d

    lp[2, 0] = offset*s3o2 # 30 degrees
    lp[2, 2] = offset*0.5

    lp[3, 0] = offset*s2o2 # 45 degrees
    lp[3, 2] = offset*s2o2

    lp[4, 0] = offset*0.5 # 60 degrees
    lp[4, 2] = offset*s3o2

    lp[5, 0] = offset*c75d # 75 degrees
    lp[5, 2] = offset*s75d

    lp[6, 0] = 0.0 # 90 degrees
    lp[6, 2] = offset

    lp[7, 0] = -offset*c75d # 105 degrees
    lp[7, 2] = offset*s75d

    lp[8, 0] = -offset*0.5 # 120 degrees
    lp[8, 2] =  offset*s3o2

    lp[9, 0] = -offset*s2o2 # 135 degrees
    lp[9, 2] =  offset*s2o2

    lp[10, 0] = -offset*s3o2 # 150 degrees
    lp[10, 2] =  offset*0.5

    lp[11, 0] = -offset*c15d # 165 degrees
    lp[11, 1] = offset*s15d

    # rotate first 12 to complete ring

    for j in range(12):
        lp[j+12, 0] = -lp[j, 0]
        lp[j+12, 2] = -lp[j, 2]

    # Repeat with +dpop offset
    for j in range(24):
        lp[j+24, :] = lp[j, :]
        lp[j+24, 5] = +dpopoffset
    
    # Repeat with -dpop offset
    for j in range(24):
        lp[j+48, :] = lp[j, :]
        lp[j+48, 5] = -dpopoffset

    bunch.checkin_particles()

    return sim

def run_element(lattice, sim, steps):

    stepper = Independent_stepper_elements(steps)
    propagator = Propagator(lattice, stepper)

    simlog = Logger(
        0, LoggerV.INFO_STEP, True
    )

    propagator.propagate(sim, simlog, 1)


# defining my own tester so I can be consistent between Synergia2
# and Synergia3
def dotest(x, y):
    if abs(x) < 1.0e-15 and abs(y) < 1.0e-15:
        return True
    else:
        if abs(x) > 1.0e-15:
            return abs(1-y/x)<1.0e-15
        else: #  abs(y) > 1.0e-15:
            return abs(1-x/y)


def prop_elem(elem_txt):
    print("reading lattice: ", elem_txt, flush=True)
    lattice = read_iota_lattice_elem(elem_txt)
    print(f'read lattice, {len(lattice.get_elements())} elements, length: {lattice.get_length()}', flush=True)

    refpart = lattice.get_reference_particle()
    print(f'reference particle, mass={refpart.get_mass()}, energy={refpart.get_total_energy()}', flush=True)

    sim1 = make_sim(refpart)
    sim4 = make_sim(refpart)

    run_element(lattice, sim1, 1)
    run_element(lattice, sim4, 5)

    bunch1 = sim1.get_bunch(0,0)
    bunch1.checkout_particles()
    bunch4 = sim4.get_bunch(0,0)
    bunch4.checkout_particles()

    assert (bunch1.get_local_num() == 24*3)
    assert (bunch4.get_local_num() == 24*3)

    lp1 = bunch1.get_particles_numpy()
    lp4 = bunch4.get_particles_numpy()

    for p in range(bunch1.get_local_num()):
        j=0; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=1; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=2; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=3; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=4; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=5; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])

def test_drift():
    prop_elem(drift_txt)

def test_m2r():
    prop_elem(m2r_txt)

def test_qa2r():
    prop_elem(qa2r_txt)

if __name__ == "__main__":
    test_drift()
    test_m2r()
    test_qa2r()
