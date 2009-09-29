#!/usr/bin/env python

import sys
sys.path.append('..')

import pyparallel_utils
import numpy
from nose.tools import *

#    const int procs = 4;
#    std::vector<int > offsets(procs), counts(procs);
#    const int num_per_proc = 5;
#    decompose_1d_raw(procs, procs * num_per_proc, offsets, counts);

def test_decompose_1d_raw():
    offsets = [0, 0, 0, 0]
    counts = [0, 0, 0, 0]
    procs = 4
    num_per_proc = 5
    a,b = pyparallel_utils.decompose_1d_raw(procs, procs * num_per_proc, offsets, counts)