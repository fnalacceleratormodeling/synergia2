import sys
sys.path.append('..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')

from mpi4py import MPI
from pyfoundation import Reference_particle, Four_momentum
from pybunch import Bunch
from pybunch import Diagnostics, Diagnostics_full2
import pyconvertors
import numpy
from nose.tools import *

mass = 2.2;
total_energy = 3.0
total_num = 100
real_num = 1.0e12
proton_charge = 1
default_s = 123.4

reference_particle = Reference_particle(mass, total_energy)
comm = MPI.COMM_WORLD
bunch = Bunch(reference_particle, proton_charge, total_num, real_num, comm)
particles = bunch.get_local_particles()
particles[:, 0:6] = numpy.random.lognormal(size=[bunch.get_total_num(), 6])
s = default_s

def test_construct():
    diagnostics = Diagnostics()

def test_construct2():
    diagnostics = Diagnostics(bunch, s)

def test_get_mean():
    diagnostics = Diagnostics(bunch, s)
    mean = diagnostics.get_mean()
    assert len(mean) == 6

def test_get_std():
    diagnostics = Diagnostics(bunch, s)
    std = diagnostics.get_mean()
    assert len(std) == 6

def test_update():
    diagnostics1 = Diagnostics()
    diagnostics1.update(bunch,s)
    diagnostics2 = Diagnostics(bunch,s)
    for i in range(0,6):
        assert diagnostics1.get_mean()[i] == diagnostics2.get_mean()[i]
        assert diagnostics1.get_std()[i] == diagnostics2.get_std()[i]

def test_construct_full2():
    diagnostics = Diagnostics_full2()

def test_construct2_full2():
    diagnostics = Diagnostics_full2(bunch, s)

def test_get_mean_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    mean = diagnostics.get_mean()
    assert len(mean) == 6

def test_get_std_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    std = diagnostics.get_mean()
    assert len(std) == 6

def test_update_full2():
    diagnostics1 = Diagnostics_full2()
    diagnostics1.update(bunch,s)
    diagnostics2 = Diagnostics_full2(bunch,s)
    for i in range(0,6):
        assert diagnostics1.get_mean()[i] == diagnostics2.get_mean()[i]
        assert diagnostics1.get_std()[i] == diagnostics2.get_std()[i]

def test_get_corr_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    corr = diagnostics.get_corr()
    assert corr.shape[0] == 6
    assert corr.shape[1] == 6

def test_get_mom2_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    mom2 = diagnostics.get_corr()
    assert mom2.shape[0] == 6
    assert mom2.shape[1] == 6

def test_get_emitx_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    emitx = diagnostics.get_emitx()

def test_get_emity_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    emity = diagnostics.get_emity()

def test_get_emitz_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    emitz = diagnostics.get_emitz()

def test_get_emitxy_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    emitxy = diagnostics.get_emitxy()

def test_get_emitxyz_full2():
    diagnostics = Diagnostics_full2(bunch, s)
    emitxyz = diagnostics.get_emitxyz()

