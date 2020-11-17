import sys
import math
sys.path.append('..')
sys.path.append('../..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')
sys.path.append('../../collective')
sys.path.append('../../simulation')

from mpi4py import MPI
from foundation import Reference_particle, Four_momentum, pconstants
from bunch import Bunch
from bunch import Diagnostics_basic, Diagnostics_full2
from utils import Commxx, Logger
from simulation import Step
from collective import Space_charge_2d_open_hockney, Space_charge_3d_open_hockney
import convertors
import numpy
from nose.tools import *

mass = 0.9;
total_energy = 5.0
total_num = 100
real_num = 1.0e10
proton_charge = 1
step_length = 1.0;

reference_particle = Reference_particle(proton_charge, mass, total_energy)
comm = Commxx()
bunch = Bunch(reference_particle, total_num, real_num, comm)
particles = bunch.get_local_particles()
particles[:, 0:6] = numpy.random.lognormal(size=[bunch.get_local_num_slots(), 6])

time_step = step_length /(reference_particle.get_beta()  * pconstants.c)
step = Step(step_length)
grid = (16, 16, 16)

def test_construct_sc2d():
    sc2d = Space_charge_2d_open_hockney(comm, grid)

def test_construct_sc3d():
    sc3d = Space_charge_3d_open_hockney(comm, grid)
    
def test_apply_sc2d():
    sc2d = Space_charge_2d_open_hockney(comm, grid)
    logger = Logger(0, False)
    sc2d.apply(bunch, time_step, step, 10, logger)

def test_apply_sc3d():
    sc3d = Space_charge_3d_open_hockney(comm, grid)
    logger = Logger(0, False)
    sc3d.apply(bunch, time_step, step, 10, logger)

