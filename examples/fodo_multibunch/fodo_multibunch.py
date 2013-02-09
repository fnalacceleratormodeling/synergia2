#!/usr/bin/env python
import sys
from math import sqrt
import synergia
from fodo_multibunch_options import opts
from mpi4py import MPI

try:
    lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")

    lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                              opts.map_order)

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

    bunches = []
    comms = synergia.utils.generate_subcomms(opts.num_bunches)
    for i in range(0, opts.num_bunches):
        bunches.append(synergia.optics.generate_matched_bunch_transverse(
                  lattice_simulator, opts.emit, opts.emit, opts.stdz, opts.dpop,
                  opts.real_particles, opts.macro_particles,
                  seed=opts.seed, comm=comms[i]))

    bunch_train = synergia.bunch.Bunch_train(bunches, opts.bunch_spacing)

    bunch_train_simulator = synergia.simulation.Bunch_train_simulator(bunch_train)

    for i in range(0, bunch_train.get_size()):
        if opts.step_tracks:
            bunch_train_simulator.add_per_step(i,
                synergia.bunch.Diagnostics_bulk_track("step_tracks_%d.h5" % i,
                                                      opts.step_tracks))
        if opts.step_full2:
            bunch_train_simulator.add_per_step(i,
                synergia.bunch.Diagnostics_full2("step_full2_%d.h5" % i))
        if opts.step_particles:
            bunch_train_simulator.add_per_step(i,
                synergia.bunch.Diagnostics_particles("step_particles_%d.h5" % i))
        if opts.turn_tracks:
            bunch_train_simulator.add_per_turn(i,
                synergia.bunch.Diagnostics_bulk_track("turn_tracks_%d.h5" % i,
                                                      opts.turn_tracks))
        if opts.turn_full2:
            bunch_train_simulator.add_per_turn(i,
                synergia.bunch.Diagnostics_full2("turn_full2_%d.h5" % i))
        if opts.turn_particles:
            bunch_train_simulator.add_per_turn(i,
                synergia.bunch.Diagnostics_particles("turn_particles_%d.h5" % i))

    propagator = synergia.simulation.Propagator(stepper)
    propagator.set_checkpoint_period(opts.checkpointperiod)
    propagator.set_concurrent_io(opts.concurrentio)
    propagator.propagate(bunch_train_simulator, opts.turns, opts.maxturns, opts.verbosity)
except Exception, e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
