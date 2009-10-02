#!/usr/bin/env python

import sys
sys.path.append('..')

from mpi4py import MPI
import pyparallel_utils
import numpy
from nose.tools import *

def test_decompose_1d_raw():
    procs = 4
    num_per_proc = 5
    offsets,counts = pyparallel_utils.decompose_1d_raw(procs, 
                                                       procs * num_per_proc)
    expected_offset = 0
    for i in range(0,procs):
        assert_equal(offsets[i],expected_offset)
        assert_equal(counts[i],num_per_proc)
        expected_offset += counts[i]

def test_decompose_1d():
    length = 23
    offsets,counts = pyparallel_utils.decompose_1d(MPI.COMM_WORLD,length)
    assert_equal(offsets[0],0)
    assert_equal(counts[0],length)

def test_decompose_1d_local():
    length = 17
    local_length = pyparallel_utils.decompose_1d_local(MPI.COMM_WORLD,length)
    assert_equal(length,local_length)