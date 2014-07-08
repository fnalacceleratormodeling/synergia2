# !/usr/bin/env synergia
import synergia

from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.utils import Commxx
from synergia.lattice import MadX_reader, Lattice
from synergia.bunch import Bunch, Diagnostics_full2
from synergia.simulation import Split_operator_stepper, \
    Split_operator_stepper_elements, Dummy_collective_operator, \
    Bunch_simulator, Propagator
from synergia.collective import Space_charge_2d_bassetti_erskine, \
    Space_charge_2d_open_hockney, \
    Space_charge_3d_open_hockney

# Get options from separate options file
from fodo_space_charge_options import opts

# Define a lattice
#     Read the lattice named "fodo" from the MadX file "fodo.madx"
lattice = synergia.lattice.MadX_reader().get_lattice("fodo", "fodo.madx")

# Define a space charge operator
if opts.space_charge == "none":
    space_charge = Dummy_collective_operator("null space charge")
elif opts.space_charge == "bassetti":
    space_charge = Space_charge_2d_bassetti_erskine()
else:
    # Set up parallel space charge solvers
    # MPI communicator, True to use simple communication avoidance
    commxx = Commxx(True)
    grid = [opts.gridx, opts.gridy, opts.gridz]
    if opts.space_charge == "2d":
        space_charge = Space_charge_2d_open_hockney(commxx, grid)
    elif opts.space_charge == "3d":
        space_charge = Space_charge_3d_open_hockney(commxx, grid)

# Define a set of simulation steps
if opts.stepper == "element":
    stepper = Split_operator_stepper_elements(lattice, opts.map_order,
                                              space_charge,
                                              opts.steps_per_element)
else:
    stepper = Split_operator_stepper(lattice, opts.map_order,
                                     space_charge, opts.num_steps)
# Define a bunch
bunch = synergia.optics.generate_matched_bunch_transverse(
              stepper.get_lattice_simulator(),
              opts.x_emit, opts.y_emit, opts.z_std, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)

# Define a bunch simulator
bunch_simulator = Bunch_simulator(bunch)

# Define a set of bunch diagnostics
#     Apply full 2nd-order moment diagnostics every step
diagnostics = Diagnostics_full2("diagnostics.h5")
bunch_simulator.add_per_step(diagnostics)

# Perform the simulation
propagator = Propagator(stepper)
propagator.propagate(bunch_simulator, opts.turns, opts.max_turns, 
                     opts.verbosity)
