#!/usr/bin/env python
import synergia
from fodo_options import opts

lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)
bunch = synergia.optics.generate_matched_bunch_transverse(
              lattice_simulator, opts.emit, opts.emit, opts.stdz, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)
stepper = synergia.simulation.Independent_stepper_elements(
                            lattice_simulator, opts.steps)
diagnostics_writer_step = synergia.bunch.Diagnostics_writer(
                                    "full2.h5",
                                    synergia.bunch.Diagnostics_full2())
diagnostics_writer_turn = synergia.bunch.Diagnostics_writer(
                                    "particles.h5",
                                    synergia.bunch.Diagnostics_particles())
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, opts.turns,
                     diagnostics_writer_step, diagnostics_writer_turn,
                     opts.verbose)
