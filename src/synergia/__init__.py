# importing mpi4py here makes HopperII happy; there doesn't seem to be
# anything wrong with doing so.
#import mpi4py
from mpi4py import MPI

#from version import __version__, version_major, version_minor, version_patch, version_tweak
#import convertors
from . import foundation
from . import utils
from . import bunch
from . import lattice
#from . import optics
from . import simulation
from . import collective

# Kokkos init
simulation.init()

# Kokkos finalize
def finalize():
    import gc
    gc.collect()
    simulation.finalize()

import atexit
atexit.register(finalize)


