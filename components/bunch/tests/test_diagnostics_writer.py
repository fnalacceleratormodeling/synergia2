import sys
sys.path.append('..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')

from mpi4py import MPI
from pyfoundation import Reference_particle, Four_momentum
from pybunch import Bunch
from pybunch import Diagnostics, Diagnostics_full2
from pybunch import Diagnostics_writer, no_diagnostics
import pyconvertors
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
comm = MPI.COMM_WORLD
bunch = Bunch(reference_particle, total_num, real_num, comm)
particles = bunch.get_local_particles()
particles[:, 0:6] = numpy.random.lognormal(size=[bunch.get_total_num(), 6])

def test_construct():
    diagnostics = Diagnostics()
    diagnostics_writer = Diagnostics_writer("test_py_construct.h5", diagnostics)

def test_construct_full2():
    diagnostics = Diagnostics_full2()
    diagnostics_writer = Diagnostics_writer("test_py_construct_full2.h5", diagnostics)

def test_construct_dummy():
    dummy = Diagnostics_writer()

def test_is_dummy():
    assert(Diagnostics_writer().is_dummy())

def test_get_diagnostics():
    diagnostics = Diagnostics()
    diagnostics_writer = Diagnostics_writer("test_py_get_diagnostics.h5", diagnostics)
    retrieved_diagnostics = diagnostics_writer.get_diagnostics()
    assert(type(retrieved_diagnostics) == type(diagnostics))

def test_write():
    diagnostics = Diagnostics(bunch)
    diagnostics_writer = Diagnostics_writer("test_py_write.h5", diagnostics)
    diagnostics_writer.write()

def test_update_and_write():
    diagnostics = Diagnostics()
    diagnostics_writer = Diagnostics_writer("test_py_update_and_write.h5", diagnostics)
    diagnostics_writer.update_and_write(bunch)
    diagnostics_writer.update_and_write(bunch)

def test_update_and_write_full2():
    diagnostics = Diagnostics_full2()
    diagnostics_writer = Diagnostics_writer("test_py_update_and_write_full2.h5", diagnostics)
    diagnostics_writer.update_and_write(bunch)
    diagnostics_writer.update_and_write(bunch)

def test_no_diagnostics():
    assert(no_diagnostics().is_dummy())
