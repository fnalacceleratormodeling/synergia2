#!/usr/bin/env python
import sys
from math import sqrt
import synergia
from fodo_options import opts
from mpi4py import MPI

try:
    lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")

    for elem in lattice.get_elements():
        if (elem.get_name() == "f"):
            elem.set_string_attribute("force_diagnostics", "true")

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

    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

    bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2("step_full2.h5"))
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("turn_full2.h5"))
    bunch_simulator.add_per_forced_diagnostics_step(
                                synergia.bunch.Diagnostics_full2("forced_full2.h5"))

    propagator = synergia.simulation.Propagator(stepper)
    propagator.set_checkpoint_period(opts.checkpointperiod)
    propagator.propagate(bunch_simulator, opts.turns, opts.maxturns, opts.verbosity)
except Exception, e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
