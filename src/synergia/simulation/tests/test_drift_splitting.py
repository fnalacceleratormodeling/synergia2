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

def read_iota_lattice_elem(e_txt):
    reader = MadX_reader()
    reader.parse(e_txt)
    lattice = reader.get_lattice('machine')

    refpart = Reference_particle(1, mp, mp+KE)
    lattice.set_reference_particle(refpart)
    return lattice


# values from the C++ test
toffs = 1.0e-3;
dpoffs = 1.0e-4;
cdtoffs =  0.05;
dpopoffs = 2.2e-4;

# create the bunch of particles matching the test particles
# in the madx run
macro_particles = 16
real_particles = 1.0e10

def make_sim(refpart):
    commxx = Commxx()
    sim = Bunch_simulator.create_single_bunch_simulator(
        refpart, macro_particles, real_particles, commxx)

    bunch = sim.get_bunch(0, 0)


    bp = bunch.get_particles_numpy()

    bp[:, 0:6] = 0
    
    # particle 0 remains at 0
    bp[:, 1] = dpoffs

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

    assert (bunch1.get_local_num() == macro_particles)
    assert (bunch4.get_local_num() == macro_particles)

    lp1 = bunch1.get_particles_numpy()
    lp4 = bunch4.get_particles_numpy()

    for p in range(1):
        j=0; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j], flush=True); assert dotest(lp1[p, j], lp4[p, j] )
        j=1; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j], flush=True); assert dotest(lp1[p, j], lp4[p, j] )
        j=2; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j], flush=True); assert dotest(lp1[p, j], lp4[p, j] )
        j=3; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j], flush=True); assert dotest(lp1[p, j], lp4[p, j] )
        j=4; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j], flush=True); assert dotest(lp1[p, j], lp4[p, j] )
        j=5; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j], flush=True); assert dotest(lp1[p, j], lp4[p, j] )

def test_drift():
    prop_elem(drift_txt)

if __name__ == "__main__":
    test_drift()
