# importing mpi4py here makes HopperII happy; there doesn't seem to be
# anything wrong with doing so.
#import mpi4py
from mpi4py import MPI

#from version import __version__, version_major, version_minor, version_patch, version_tweak
#import convertors
import foundation
import utils
import bunch
import lattice
#import optics
import simulation
import collective

# Kokkos init
simulation.init()

# Kokkos finalize
def finalize():
    import gc
    gc.collect()
    simulation.finalize()

import atexit
atexit.register(finalize)


