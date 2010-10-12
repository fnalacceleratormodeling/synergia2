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
from pylattice import Lattice, xml_save_lattice, xml_load_lattice
from pysimulation import Collective_operator, Lattice_simulator, \
    Split_operator_stepper, Propagator
from pybunch import Diagnostics_full2, Diagnostics_particles, \
    Diagnostics_writer, no_diagnostics
from matching import get_matched_bunch_transverse_parameters
from pyconvertors import xml_save_array1d, xml_save_array2d
from one_turn_map import linear_one_turn_map
import numpy
import sys

num_macro_particles = 32000
seed = 4
grid = [16, 16, 16]
num_real_particles = 1e12
num_steps = 8
num_turns = 4
map_order = 2
emit = 1e-6
stdz = 0.01
dpop = 1e-4

lattice = Mad8_reader().get_lattice("fodo", "fodo.lat")
xml_save_lattice(lattice, "cxx_lattice.xml")
lattice_simulator = Lattice_simulator(lattice, map_order)
means, covariance_matrix = \
    get_matched_bunch_transverse_parameters(lattice_simulator,
                                            emit, emit, stdz, dpop)
xml_save_array1d(means,"cxx_means.xml")
xml_save_array2d(covariance_matrix,"cxx_covariance_matrix.xml")
