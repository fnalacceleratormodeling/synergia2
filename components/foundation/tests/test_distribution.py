#!/usr/bin/env python

import sys
sys.path.append('..')

from mpi4py import MPI
from pyfoundation import Distribution, Random_distribution
import numpy
from nose.tools import *

test_seed = 12345678
    
def test_construct():
    r = Random_distribution(test_seed,MPI.COMM_WORLD)
