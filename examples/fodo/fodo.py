#!/usr/bin/env python
import sys
import synergia
from fodo_options import opts

lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")

# Set the same aperture radius for all elements
for elem in lattice.get_elements():
    elem.set_double_attribute("aperture_radius", opts.radius)

lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                          opts.map_order)

bunch = synergia.optics.generate_matched_bunch_transverse(
              lattice_simulator, opts.emit, opts.emit, opts.stdz, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)

if opts.stepper == "splitoperator":
    # Use the Split operator stepper with a dummy collective operator
    # (with evenly-spaced steps)
    no_op = synergia.simulation.Dummy_collective_operator("stub")
    stepper = synergia.simulation.Split_operator_stepper(
                            lattice_simulator, no_op, opts.steps)
elif opts.stepper == "independent":
    # Use the Independent particle stepper (by element)
    stepper = synergia.simulation.Independent_stepper_elements(
                            lattice_simulator, opts.steps)
else:
    sys.stderr.write("fodo.py: stepper must be either 'independent' or 'splitoperator'\n")
    sys.exit(1)

multi_diagnostics_writer_step = synergia.bunch.Multi_diagnostics_writer()
for part in range(0, opts.step_tracks):
    multi_diagnostics_writer_step.append(synergia.bunch.Diagnostics_writer(
                                    "step_track_%02d.h5" % part,
                                    synergia.bunch.Diagnostics_track(part)))
if opts.step_full2:
    multi_diagnostics_writer_step.append(synergia.bunch.Diagnostics_writer(
                                    "step_full2.h5",
                                    synergia.bunch.Diagnostics_full2()))
if opts.step_particles:
    multi_diagnostics_writer_step.append(synergia.bunch.Diagnostics_writer(
                                    "step_particles.h5",
                                    synergia.bunch.Diagnostics_particles()))

multi_diagnostics_writer_turn = synergia.bunch.Multi_diagnostics_writer()
for part in range(0, opts.turn_tracks):
    multi_diagnostics_writer_turn.append(synergia.bunch.Diagnostics_writer(
                                    "turn_track_%02d.h5" % part,
                                    synergia.bunch.Diagnostics_track(part)))
if opts.turn_full2:
    multi_diagnostics_writer_turn.append(synergia.bunch.Diagnostics_writer(
                                    "turn_full2.h5",
                                    synergia.bunch.Diagnostics_full2()))
if opts.turn_particles:
    multi_diagnostics_writer_turn.append(synergia.bunch.Diagnostics_writer(
                                    "turn_particles.h5",
                                    synergia.bunch.Diagnostics_particles()))

propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, opts.turns,
                     multi_diagnostics_writer_step,
                     multi_diagnostics_writer_turn,
                     opts.verbose)
