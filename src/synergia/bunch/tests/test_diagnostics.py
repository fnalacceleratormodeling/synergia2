import sys
sys.path.append('..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')

from mpi4py import MPI
from foundation import Reference_particle, Four_momentum
from bunch import Bunch
from bunch import Diagnostics_basic, Diagnostics_full2
import convertors
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
    diagnostics = Diagnostics_basic(bunch, "dummy.h5")

def test_get_s():
    diagnostics = Diagnostics_basic(bunch, "dummy.h5")
    diagnostics.update()
    assert_almost_equal(diagnostics.get_s(), partial_s)

def test_get_repetition():
    diagnostics = Diagnostics_basic(bunch, "dummy.h5")
    diagnostics.update()
    assert_equal(diagnostics.get_repetition(), turns)

def test_get_trajectory_length():
    diagnostics = Diagnostics_basic(bunch, "dummy.h5")
    diagnostics.update()
    assert_almost_equal(diagnostics.get_trajectory_length(),
                        turns * turn_length + partial_s)

def test_get_mean():
    diagnostics = Diagnostics_basic(bunch, "dummy.h5")
    diagnostics.update()
    mean = diagnostics.get_mean()
    assert len(mean) == 6

def test_get_std():
    diagnostics = Diagnostics_basic(bunch, "dummy.h5")
    diagnostics.update()
    std = diagnostics.get_mean()
    assert len(std) == 6

def test_update():
    diagnostics1 = Diagnostics_basic(bunch, "dummy1.h5")
    diagnostics1.update()
    diagnostics2 = Diagnostics_basic(bunch, "dummy2.h5")
    diagnostics2.update()
    for i in range(0, 6):
        assert diagnostics1.get_mean()[i] == diagnostics2.get_mean()[i]
        assert diagnostics1.get_std()[i] == diagnostics2.get_std()[i]

def test_construct_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")

def test_get_s_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    assert_almost_equal(diagnostics.get_s(), partial_s)

def test_get_repetition_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    assert_equal(diagnostics.get_repetition(), turns)

def test_get_trajectory_length_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    assert_almost_equal(diagnostics.get_trajectory_length(),
                        turns * turn_length + partial_s)

def test_get_mean_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    mean = diagnostics.get_mean()
    assert len(mean) == 6

def test_get_std_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    std = diagnostics.get_mean()
    assert len(std) == 6

def test_update_full2():
    diagnostics1 = Diagnostics_full2(bunch, "dummy1.h5")
    diagnostics1.update()
    diagnostics2 = Diagnostics_full2(bunch, "dummy2.h5")
    diagnostics2.update()
    for i in range(0, 6):
        assert diagnostics1.get_mean()[i] == diagnostics2.get_mean()[i]
        assert diagnostics1.get_std()[i] == diagnostics2.get_std()[i]

def test_get_corr_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    corr = diagnostics.get_corr()
    assert corr.shape[0] == 6
    assert corr.shape[1] == 6

def test_get_mom2_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    mom2 = diagnostics.get_corr()
    assert mom2.shape[0] == 6
    assert mom2.shape[1] == 6

def test_get_emitx_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    emitx = diagnostics.get_emitx()

def test_get_emity_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    emity = diagnostics.get_emity()

def test_get_emitz_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    emitz = diagnostics.get_emitz()

def test_get_emitxy_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    emitxy = diagnostics.get_emitxy()

def test_get_emitxyz_full2():
    diagnostics = Diagnostics_full2(bunch, "dummy.h5")
    diagnostics.update()
    emitxyz = diagnostics.get_emitxyz()

