# !/usr/bin/env synergia
import synergia

from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Mad8_reader, Lattice
from synergia.bunch import Bunch, Diagnostics_basic 
from synergia.simulation import Independent_stepper_elements, Bunch_simulator, \
    Propagator
from fodo_options import opts
from mpi4py import MPI
import sys

# We wrap the entire simulation in a try..except block in order to allow
# for graceful failures under MPI.
try:
    # Read the lattice named "fodo" from the Mad8 file "fodo.lat"
    lattice = Mad8_reader().get_lattice("fodo", "fodo.lat")
    
    # put an aperture on the defocussing quadrupole to generate a loss
    for elem in lattice.get_elements():
        if elem.get_name() == "d":
            elem.set_string_attribute("aperture_type", "rectangular")
            # vertical aperture is .002 above and below
            elem.set_double_attribute("rectangular_aperture_height", 2*0.002)
            elem.set_double_attribute("rectangular_aperture_width", 10.0)

    # Define the simulation steps
    stepper = Independent_stepper_elements(lattice, opts.map_order, 
                                           opts.steps_per_element)
    
    # Define the parameters for the bunch
    bunch = synergia.optics.generate_matched_bunch_transverse(
                  stepper.get_lattice_simulator(),
                  opts.x_emit, opts.y_emit, opts.z_std, opts.dpop,
                  opts.real_particles, opts.macro_particles,
                  seed=opts.seed)
    
    bunch.set_bucket_index(0)

    # create the loss diagnostics
    loss_diag = synergia.lattice.Diagnostics_loss("loss_diagnostics.h5", "aperture")
    # set the bunch for the loss
    loss_diag.set_bunch(bunch)

    # add the loss diagnostics to the lattice
    lattice.add_loss_diagnostics(loss_diag)

    # Apply basic diagnostics every step
    diagnostics = Diagnostics_basic("diagnostics.h5")
    bunch_simulator = Bunch_simulator(bunch)
    bunch_simulator.add_per_step(diagnostics)
    
    # Perform the simulation
    propagator = Propagator(stepper)
    propagator.propagate(bunch_simulator, opts.turns, opts.max_turns, 
                         opts.verbosity)

except Exception, e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
