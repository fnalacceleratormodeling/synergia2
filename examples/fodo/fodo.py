#!/usr/bin/env python
import sys
from math import sqrt
import synergia
from fodo_options import opts
from mpi4py import MPI

try:
    lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")

    # Set the same aperture radius for all elements
    for elem in lattice.get_elements():
        elem.set_double_attribute("circular_aperture_radius", opts.radius)

    lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                              opts.map_order)

    if opts.elliptical:
        aperture_operation_extractor_map = lattice_simulator.get_aperture_operation_extractor_map()
        aperture_operation_extractor_map.set_extractor("default",
                                                       synergia.simulation.Elliptical_extractor())
        lattice.set_all_double_attribute("elliptical_aperture_horizontal_radius",
                                     opts.radius*sqrt(opts.aspect))
        lattice.set_all_double_attribute("elliptical_aperture_vertical_radius",
                                     opts.radius/sqrt(opts.aspect))

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

    if opts.step_tracks:
        bunch_simulator.add_per_step(synergia.bunch.Diagnostics_bulk_track("step_tracks.h5", opts.turn_tracks/MPI.COMM_WORLD.size))

    if opts.step_full2:
        bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2("step_full2.h5"))
    if opts.step_particles:
        bunch_simulator.add_per_step(synergia.bunch.Diagnostics_particles("step_particles.h5"))

    if opts.turn_tracks:
        bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("turn_tracks.h5", opts.turn_tracks/MPI.COMM_WORLD.size))

    if opts.turn_full2:
        bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("turn_full2.h5"))
    if opts.turn_particles:
        bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("turn_particles.h5"))

    propagator = synergia.simulation.Propagator(stepper)
    propagator.set_checkpoint_period(opts.checkpointperiod)
    propagator.set_concurrent_io(opts.concurrentio)
    propagator.propagate(bunch_simulator, opts.turns, opts.maxturns, opts.verbosity)
except Exception, e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
