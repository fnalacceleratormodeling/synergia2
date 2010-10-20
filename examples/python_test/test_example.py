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
space_charge = synergia.simulation.Collective_operator("space charge")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, map_order)
stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, space_charge,
                                          num_steps)
bunch = synergia.optics.generate_matched_bunch_transverse(lattice_simulator, emit, emit, stdz, dpop,
                                        num_real_particles, num_macro_particles,
                                        seed=seed)
diagnostics_writer_step = synergia.bunch.Diagnostics_writer("example_full2.h5",
                                                            synergia.bunch.Diagnostics_full2())
diagnostics_writer_turn = synergia.bunch.Diagnostics_writer("example_particles.h5",
                                                            synergia.bunch.Diagnostics_particles())
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, num_turns, diagnostics_writer_step, diagnostics_writer_turn)
