# !/usr/bin/env synergia
import synergia

from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Mad8_reader, Lattice
from synergia.bunch import Bunch, Diagnostics_basic, Diagnostics_track
from synergia.simulation import Independent_stepper_elements, Bunch_simulator, \
    Propagator

# Get options from separate options file
from fodo_workflow_options import opts

# Define a lattice
#     Read the lattice named "fodo" from the Mad8 file "fodo.lat"
lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")

if opts.libff:
    extractor_type = "libff"
else:
    extractor_type = "chef_propagate"
for element in lattice.get_elements():
    element.set_string_attribute("extractor_type", extractor_type)

# Define a set of simulation steps
stepper = Independent_stepper_elements(lattice, opts.map_order, 
                                       opts.steps_per_element)

# Define a bunch
bunch = synergia.optics.generate_matched_bunch_transverse(
              stepper.get_lattice_simulator(),
              opts.x_emit, opts.y_emit, opts.z_std, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)

# Define a bunch simulator
bunch_simulator = Bunch_simulator(bunch)

# Define a set of bunch diagnostics
#     Apply basic diagnostics every step
diagnostics = Diagnostics_basic("diagnostics.h5")
bunch_simulator.add_per_step(diagnostics)

# Optionally, save track of a single particle
if opts.track_particle:
    track_diagnostics = Diagnostics_track("track.h5", opts.track_particle)
    bunch_simulator.add_per_step(track_diagnostics)

# Perform the simulation
propagator = Propagator(stepper)
propagator.propagate(bunch_simulator, opts.turns, opts.max_turns, 
                     opts.verbosity)
