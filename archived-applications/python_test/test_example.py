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
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, map_order)
bunch = synergia.optics.generate_matched_bunch_transverse(lattice_simulator, emit, emit, stdz, dpop,
                                        num_real_particles, num_macro_particles,
                                        seed=seed)
space_charge = synergia.collective.Space_charge_3d_open_hockney(bunch.get_comm(), grid)
stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, space_charge,
                                          num_steps)
diagnostics_step = synergia.bunch.Diagnostics_full2(bunch, "example_full2.h5")
diagnostics_turn = synergia.bunch.Diagnostics_particles(bunch, "example_particles.h5")
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, num_turns, diagnostics_step, diagnostics_turn)
