#!/usr/bin/env python

import sys
sys.path.append('..')

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
        assert_equal(counts[i],num_per_proc)
        assert_equal(offsets[i],expected_offset)
        expected_offset += counts[i]
