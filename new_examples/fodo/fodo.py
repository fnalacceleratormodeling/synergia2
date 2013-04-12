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
    
    # Define the simulation steps
    stepper = Independent_stepper_elements(lattice, opts.map_order, 
                                           opts.steps_per_element)
    
    # Define the parameters for the bunch
    bunch = synergia.optics.generate_matched_bunch_transverse(
                  stepper.get_lattice_simulator(),
                  opts.x_emit, opts.y_emit, opts.z_std, opts.dpop,
                  opts.real_particles, opts.macro_particles,
                  seed=opts.seed)
    
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
