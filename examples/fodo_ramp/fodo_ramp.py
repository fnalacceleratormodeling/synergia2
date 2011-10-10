#!/usr/bin/env python
import sys
import synergia
from fodo_ramp_options import opts

class Ramp_actions(synergia.simulation.Propagate_actions):
    def __init__(self, multiplier):
        synergia.simulation.Propagate_actions.__init__(self)
        self.multiplier = multiplier
    def turn_end_action(self, stepper, bunch, turn_num):
        print "modifying lattice"
        for element in stepper.get_lattice_simulator().get_lattice().get_elements():
            if element.get_type() == "quadrupole":
                old_k1 = element.get_double_attribute("k1")
                element.set_double_attribute("k1", self.multiplier*old_k1)
                print "  updated", element.get_name()
        stepper.get_lattice_simulator().update()


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

diagnostics_actions = synergia.simulation.Standard_diagnostics_actions()
for part in range(0, opts.step_tracks):
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_track(bunch,
                                                                   "step_track_%02d.h5" % part,
                                                                   part))
if opts.step_full2:
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_full2(bunch, "step_full2.h5"))
if opts.step_particles:
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_particles(bunch,
                                                                       "step_particles.h5"))

for part in range(0, opts.turn_tracks):
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_track(bunch,
                                                                   "turn_track_%02d.h5" % part,
                                                                   part))
if opts.turn_full2:
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_full2(bunch, "turn_full2.h5"))
if opts.turn_particles:
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_particles(
                                    bunch, "turn_particles.h5"))

ramp_actions = Ramp_actions(1.1)

propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, opts.turns,
                     diagnostics_actions,
                     ramp_actions,
                     opts.verbosity)
