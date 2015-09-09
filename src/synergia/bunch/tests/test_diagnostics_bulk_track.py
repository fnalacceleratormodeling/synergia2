import sys
sys.path.append('..')
sys.path.append('../..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')
sys.path.append('../../lattice')
sys.path.append('../../simulation')

from mpi4py import MPI
from foundation import Reference_particle, Four_momentum
from bunch import Bunch
from bunch import Diagnostics_bulk_track
from utils import Commxx
import convertors
import lattice
import simulation
import numpy
from nose.tools import *

mass = 2.2;
total_energy = 3.0
total_num = 100
real_num = 1.0e12
proton_charge = 1
turns = 17;
turn_length = 246.8;
partial_s = 123.4;


reference_particle = Reference_particle(proton_charge, mass, total_energy)
reference_particle.set_trajectory(turns, turn_length, partial_s)
comm = Commxx()
bunch = Bunch(reference_particle, total_num, real_num, comm)
particles = bunch.get_local_particles()
particles[:, 0:6] = numpy.random.lognormal(size=[bunch.get_total_num(), 6])

def test_construct():
    diagnostics = Diagnostics_bulk_track("dummy.h5", 0)

def test_construct2():
    diagnostics = Diagnostics_bulk_track("dummy.h5", 0)

def test_update():
    diagnostics = Diagnostics_bulk_track("dummy.h5", 0)
    diagnostics.set_bunch(bunch)
    diagnostics.update()
