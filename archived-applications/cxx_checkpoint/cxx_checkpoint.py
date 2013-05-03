#!/usr/bin/env python
import synergia

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

lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")
synergia.lattice.xml_save_lattice(lattice, "cxx_lattice.xml")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, map_order)
means, covariance_matrix = \
    synergia.optics.get_matched_bunch_transverse_parameters(lattice_simulator,
                                            emit, emit, stdz, dpop)
synergia.convertors.xml_save_array1d(means,"cxx_means.xml")
synergia.convertors.xml_save_array2d(covariance_matrix,"cxx_covariance_matrix.xml")
