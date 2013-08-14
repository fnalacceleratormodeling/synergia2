# !/usr/bin/env synergia
import synergia

from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Mad8_reader, Lattice
from synergia.bunch import Bunch, Diagnostics_basic 
from synergia.simulation import Independent_stepper_elements, Bunch_simulator, \
    Propagator

# Define a lattice
#     Read the lattice named "fodo" from the Mad8 file "fodo.lat"
lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")


# Define a set of simulation steps
map_order = 1
steps_per_element = 2
stepper = Independent_stepper_elements(lattice, map_order, steps_per_element)

# Define a bunch
x_emit = 1.0e-6  # m-rad, RMS
y_emit = 1.0e-6  # m-rad, RMS
z_std = 0.01  # m
dpop = 1.0e-4  # unitless, RMS \frac{\delta p}{p_{tot}}
real_particles = 1.2e12  # unitless, meaningless in this simulation
                         #           without collective effects
macro_particles = 100
seed = 1415926  # random number seed; 0 for automatic calculation (GSL)
bunch = synergia.optics.generate_matched_bunch_transverse(
              stepper.get_lattice_simulator(),
              x_emit, y_emit, z_std, dpop,
              real_particles, macro_particles,
              seed=seed)

# Define a bunch simulator
bunch_simulator = Bunch_simulator(bunch)

# Define a set of bunch diagnostics
#     Apply basic diagnostics every step
diagnostics = Diagnostics_basic("diagnostics.h5")
bunch_simulator.add_per_step(diagnostics)

# Perform the simulation
propagator = Propagator(stepper)
turns = 4  # really repetitions, since this isn't a ring
max_turns = 0 # Number of turns to run before writing checkpoint and stopping
              # When max_turns is 0, the simulation continues until the end.
verbosity = 2  # Display information about each simulation step
propagator.propagate(bunch_simulator, turns, max_turns, verbosity)
