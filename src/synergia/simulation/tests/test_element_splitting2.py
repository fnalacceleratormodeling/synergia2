#!/usr/bin/env python

import sys
#sys.path.append('../../..')
#import local_paths


from nose.tools import *
from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice, MadX_reader
from synergia.bunch import Bunch
from synergia.utils import Commxx
from synergia.utils.parallel_utils import Logger
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
    lattice.set_all_string_attribute('extractor_type', 'libff')

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

def make_bunch(refpart):
    commxx = Commxx()
    bunch = Bunch(refpart, macro_particles, real_particles, commxx)


    bp = bunch.get_local_particles()

    bp[:, 0:6] = 0
    
    # particle 0 remains at 0
    bp[1, 0] = toffs # particle 1 offset in x
    bp[2, 1] = dpoffs # particle 2 momentum offset
    bp[3, 2] = toffs
    bp[4, 3] = dpoffs
    bp[5, 4] = cdtoffs
    bp[6, 5] = dpopoffs

    bp[7, 0] = -toffs
    bp[7, 5] = dpopoffs

    bp[8, 1] = -dpoffs
    bp[8, 5] = dpopoffs

    bp[9, 2] = toffs
    bp[9, 3] = -dpoffs
    bp[9, 5] = -dpopoffs

    bp[10, 0] = -toffs
    bp[10, 1] = -dpoffs
    bp[10, 2] = toffs
    bp[10, 3] = dpoffs
    bp[10, 5] = dpopoffs

    bp[11, 0] = toffs
    bp[11, 1] = dpoffs
    bp[11, 2] = -toffs
    bp[11, 3] = -dpoffs
    bp[11, 5] = -dpopoffs

    return bunch

def run_element(lattice, bunch, steps):

    sim = Bunch_simulator(bunch)
    stepper = Independent_stepper_elements(lattice, 1, steps)
    propagator = Propagator(stepper)

    propagator.propagate(sim, 1)

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
    #lattice.set_all_string_attribute('extractor_type', 'chef_propagate')
    lattice.set_all_string_attribute('extractor_type', 'libff')
    print(f'read lattice, {len(lattice.get_elements())} elements, length: {lattice.get_length()}', flush=True)

    refpart = lattice.get_reference_particle()
    print(f'reference particle, mass={refpart.get_mass()}, energy={refpart.get_total_energy()}', flush=True)

    bunch1 = make_bunch(refpart)
    bunch4 = make_bunch(refpart)

    run_element(lattice, bunch1, 1)
    run_element(lattice, bunch4, 5)

    assert (bunch1.get_local_num() == macro_particles)
    assert (bunch4.get_local_num() == macro_particles)

    lp1 = bunch1.get_local_particles()
    lp4 = bunch4.get_local_particles()

    for p in range(bunch1.get_local_num()):
        j=0; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=1; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=2; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=3; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=4; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])
        j=5; print("particle ", p, " coordinate ", j, lp1[p,j], ' <--> ', lp4[p,j]); assert dotest( lp1[p, j], lp4[p, j])

        #  uncomment the next two lines to provoke a failure
        # if p == bunch1.get_local_num()-1 and j == 5:
        # lp1[p, j] += lp1[p, j] + 1
        # assert_almost_equal( lp1[p, j], lp4[p, j], 13)

def test_drift():
    prop_elem(drift_txt)

def test_m1r():
    prop_elem(m1r_txt)

def test_m2r():
    prop_elem(m2r_txt)

def test_qa2r():
    prop_elem(qa2r_txt)

if __name__ == "__main__":
    test_drift()
    test_m1r()
    test_m2r()
    test_qa2r()

