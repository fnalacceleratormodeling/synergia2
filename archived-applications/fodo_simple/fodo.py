#!/usr/bin/env python
import synergia
from fodo_options import opts

lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)
bunch = synergia.optics.generate_matched_bunch_transverse(
              lattice_simulator, opts.emit, opts.emit, opts.stdz, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)


diagnostics_step = synergia.bunch.Diagnostics_full2("full2.h5")
diagnostics_turn = synergia.bunch.Diagnostics_particles("particles.h5")

bunch_simulator.add_per_step(diagnostics_step)
bunch_simulator.add_per_turn(diagnostics_turn)

stepper = synergia.simulation.Independent_stepper_elements(
                            lattice_simulator, opts.steps)


propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch_simulator, opts.turns, opts.turns, opts.verbose)
