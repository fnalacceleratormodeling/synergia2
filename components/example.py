#!/usr/bin/env python

# jfa; This file is a goal version of a simple script. It does not (yet) work,
# and never may.

import synergia

num_macro_particles = 1000
seed = 4
grid = [16, 16, 16]
num_real_particles = 1e12
num_steps = 10
num_turns = 1000

lattice = synergia.Mad8_reader().read("fodo.lat")
bunch = synergia.generate_matched_bunch(lattice, num_particles, num_real_particles, seed=seed)
space_charge = synergia.Space_charge_3d_open_hockney(grid)
stepper = synergia.Split_operator_stepper(lattice, num_steps)
propagator = synergia.Propagator(stepper)
propagator.propagate(bunch, num_turns, diagnostics_per_step=False, diagnostics_per_turn=True)
