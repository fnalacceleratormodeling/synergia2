#!/usr/bin/env python

# script to perform same operations in python as cxx_benchmark does in
# C++

import sys, os
import numpy as np
import synergia
import mpi4py.MPI

from benchmark_options import opts

grid_shape = [opts.gridx, opts.gridy, opts.gridz]
part_per_cell = opts.partpercell
num_macro_particles = grid_shape[0] * grid_shape[1] * grid_shape[2] * part_per_cell

seed = 4
num_real_particles = 1.0e13
num_steps = 8
num_turns = 4
map_order = 2

lattice = synergia.lattice.Lattice()
synergia.lattice.xml_load_lattice(lattice, "cxx_lattice.xml")
lattice.set_all_string_attribute("extractor_type", "chef_propagate", False)

comm = synergia.utils.Commxx()
bunch = synergia.bunch.Bunch(lattice.get_reference_particle(), num_macro_particles, num_real_particles, comm)

distribution = synergia.foundation.Random_distribution(seed, comm)

means = np.zeros((6), dtype='d')
covariances = np.zeros((6,6), dtype='d')

means = np.loadtxt("cxx_means.txt")
covariances = np.loadtxt("cxx_covariance_matrix.txt")

synergia.bunch.populate_6d(distribution, bunch, means, covariances)
if opts.sortperiod > 0:
    bunch.sort(4)
bunch.set_sort_period(opts.sortperiod)

space_charge = synergia.collective.Space_charge_3d_open_hockney(grid_shape)

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, map_order)

stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, space_charge, num_steps)

propagator = synergia.simulation.Propagator(stepper)

propagator.set_final_checkpoint(False)

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

if opts.diagnostics:
    bunch_simulator.add_per_step(synergia.bunch.Diagnostics_basic("diag_per_step.h5"))
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("diag_per_turn.h5"))

max_turns = 0
t0 = mpi4py.MPI.Wtime()
propagator.propagate(bunch_simulator, num_turns, max_turns, opts.verbosity)
t1 = mpi4py.MPI.Wtime()
print "propagate time = ", t1-t0
