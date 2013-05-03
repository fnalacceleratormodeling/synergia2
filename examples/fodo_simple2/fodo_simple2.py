# !/usr/bin/env synergia
import synergia

from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Lattice_element, Lattice, Mad8_adaptor_map
from synergia.bunch import Bunch, Diagnostics_basic 
from synergia.simulation import Independent_stepper_elements, Bunch_simulator, \
    Propagator

# Set the lattice element parameters
focus_length = 7  # meters
quad_sep = 10  # meters
quad_length = 2.0  # meters
strength = 1.0 / (focus_length * quad_length)  # 1/meters^2

# Define the lattice elements
o = Lattice_element("drift", "o")
o.set_double_attribute("l", quad_sep - quad_length)
f = Lattice_element("quadrupole", "f")
f.set_double_attribute("l", quad_length)
f.set_double_attribute("k1", strength)
d = Lattice_element("quadrupole", "d")
d.set_double_attribute("l", quad_length)
d.set_double_attribute("k1", -strength)

# Define the  fodo lattice itself, interpreting elements
# by their Mad8 definitions
lattice = Lattice("fodo", Mad8_adaptor_map())
# Add copies of the lattice elements to the fodo lattice
lattice.append(f)
lattice.append(o)
lattice.append(d)
lattice.append(o)

# Define a reference particle
total_energy = 1.5  # GeV
four_momentum = Four_momentum(pconstants.proton_mass, total_energy)
reference_particle = Reference_particle(pconstants.proton_charge,
                                        four_momentum)
lattice.set_reference_particle(reference_particle)

# Define the simulation steps
map_order = 1
steps_per_element = 2
stepper = Independent_stepper_elements(lattice, map_order, steps_per_element)

# Define the parameters for the bunch
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

# Apply basic diagnostics every step
diagnostics = Diagnostics_basic("diagnostics.h5")
bunch_simulator = Bunch_simulator(bunch)
bunch_simulator.add_per_step(diagnostics)

# Perform the simulation
propagator = Propagator(stepper)
turns = 4  # really repetitions, since this isn't a ring
max_turns = 0 # Number of turns to run before writing checkpoint and stopping
              # When max_turns is 0, the simulation continues until the end.
verbosity = 2  # Display information about each simulation step
propagator.propagate(bunch_simulator, turns, max_turns, verbosity)
