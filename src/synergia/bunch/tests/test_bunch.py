#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.append('..')
sys.path.append('../..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')

from mpi4py import MPI
from foundation import Reference_particle, Four_momentum
from bunch import Bunch, Fixed_t_z_zeroth
from utils import Commxx
import convertors
import numpy
from nose.tools import *

mass = 2.2;
total_energy = 3.0
total_num = 100
real_num = 1.0e12
proton_charge = 1
reference_particle = Reference_particle(proton_charge, mass, total_energy)

def test_construct():
    b = Bunch(reference_particle, total_num, real_num,
               Commxx())

def test_construct_periodic():
    b = Bunch(reference_particle, total_num, real_num, Commxx())
    z_period_length=4.6
    b.set_z_period_length(z_period_length)
    assert_equal(b.get_z_period_length(),z_period_length)
    assert_equal(b.is_periodic(),True)

def test_construct_longitudinal_aperture():
    b = Bunch(reference_particle, total_num, real_num, Commxx())
    z_length=4.6
    b.set_longitudinal_aperture_length(z_length)
    assert_equal(b.get_longitudinal_aperture_length(),z_length)
    assert_equal(b.has_longitudinal_aperture(),True)

def test_construct_bucket_index():
    b = Bunch(reference_particle, total_num, real_num, Commxx())
    assert_equal(b.is_bucket_index_assigned(),False)
    bk=2
    b.set_bucket_index(bk)
    assert_equal(b.get_bucket_index(),bk)
    assert_equal(b.is_bucket_index_assigned(),True)

def test_set_get_particle_charge():
    b = Bunch(reference_particle, total_num, real_num,
               Commxx())
    new_charge = -proton_charge
    b.set_particle_charge(new_charge)
    assert_equal(b.get_particle_charge(),new_charge)

def test_set_get_real_num():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    new_real_num = 2*real_num
    b.set_real_num(new_real_num)
    assert_equal(b.get_real_num(),new_real_num)

def test_set_get_local_num():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    new_local_num = b.get_local_num() - 3
    b.set_local_num(new_local_num)
    assert_equal(b.get_local_num(),new_local_num)

def test_update_total():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    new_local_num = b.get_local_num() - 5
    new_total_num = Commxx().get_size() * new_local_num
    b.set_local_num(new_local_num)
    b.update_total_num()
    assert_equal(b.get_total_num(),new_total_num)

def test_convert_to_state():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    b.convert_to_state(Bunch.fixed_t_lab)


def test_get_reference_particle():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    r = b.get_reference_particle()
    assert_equal(reference_particle.get_total_energy(),r.get_total_energy())

def test_get_local_particles():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    p = b.get_local_particles()
    assert p.shape[0] >= b.get_local_num()
    assert p.shape[1] == 7

def test_get_mass():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    assert_equal(reference_particle.get_mass(),
                 b.get_mass())

def test_get_state():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())
    assert b.get_state() == Bunch.fixed_z_lab

def test_get_comm():
    b = Bunch(reference_particle, total_num, real_num,
              Commxx())

def test_attrs():
    assert Bunch.x == 0
    assert Bunch.xp == 1
    assert Bunch.y == 2
    assert Bunch.yp == 3
    assert Bunch.z == 4
    assert Bunch.zp == 5
    assert Bunch.cdt == 4
    assert Bunch.dpop == 5
    assert Bunch.id == 6

