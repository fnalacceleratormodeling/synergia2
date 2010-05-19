#!/usr/bin/env python


import sys
sys.path.append('foundation')
sys.path.append('lattice')
sys.path.append('simulation')
sys.path.append('bunch')
sys.path.append('optics')
sys.path.append('convertors')
sys.path.append("/home/amundson/work/synergia2-devel_1_0/install/lib")

from mad8_reader import Mad8_reader
from pysimulation import Collective_operator, Lattice_simulator, \
    Split_operator_stepper, Propagator
from pybunch import Diagnostics_full2, Diagnostics_writer, no_diagnostics
from matching import generate_matched_bunch_transverse
import pyconvertors

num_macro_particles = 1000
seed = 4
grid = [16, 16, 16]
num_real_particles = 1e12
num_steps = 4
num_turns = 2
map_order = 2
emit = 1e-6
stdz = 0.01
dpop = 1e-4

lattice = Mad8_reader().get_lattice("fodo", "optics/tests/fodo.lat")
#space_charge = Space_charge_3d_open_hockney(grid)
space_charge = Collective_operator("space charge")
lattice_simulator = Lattice_simulator(lattice, map_order)
stepper = Split_operator_stepper(lattice_simulator, space_charge,
                                          num_steps)
bunch = generate_matched_bunch_transverse(lattice_simulator, emit, emit, stdz, dpop,
                                        num_real_particles, num_macro_particles,
                                        seed=seed)
diagnostics = Diagnostics_full2()
diagnostics_writer = Diagnostics_writer("example_full2.h5", diagnostics)
propagator = Propagator(stepper)
propagator.propagate(bunch, num_turns, diagnostics_writer, no_diagnostics())
