#!/usr/bin/env python

import sys
sys.path.append('..')
sys.path.append('../../convertors')

from mpi4py import MPI
from foundation import Distribution, Random_distribution
import convertors
import numpy
from nose.tools import *

test_seed = 12345678
array_length = 1000

def test_construct():
    r = Random_distribution(0,MPI.COMM_WORLD)

def test_construct2():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)

def test_construct3():
    r = Random_distribution(test_seed,MPI.COMM_WORLD,Random_distribution.mt19937)

def test_get_seed():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)
    assert_equal(test_seed,r.get_original_seed())

def test_set_seed():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)
    new_seed = 987654321
    r.set_seed(new_seed)
    assert_equal(new_seed,r.get_original_seed())

def test_get_default_seed():
    assert(Random_distribution.get_default_seed("/dev/urandom") != 0)

def test_get():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)
    assert(r.get() != r.get())

def test_fill_uniform():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)
    array = numpy.zeros([array_length],'d')
    min = 3.4
    max = 7.79
    r.fill_uniform(array,min,max)

def test_fill_unit_gaussian():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)
    array = numpy.zeros([array_length],'d')
    r.fill_unit_gaussian(array)

def test_fill_unit_disk():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)
    x = numpy.zeros([array_length],'d')
    y = numpy.zeros([array_length],'d')
    r.fill_unit_disk(x,y)

