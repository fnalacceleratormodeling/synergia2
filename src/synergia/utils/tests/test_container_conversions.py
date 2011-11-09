#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from mpi4py import MPI
import synergia
import numpy
from nose.tools import *
from wrap_containers_python  import *
from convertors import *

def list_to_std_vector_int():
    pylist=[4,3,6]
    take_vector_int(pylist)

def list_to_std_list_int():
    pylist=[14,13,16]
    take_list_int(pylist)
    
def list_to_std_vector_dd():
    pylist=[4.,3.,6.]
    take_vector_dd(pylist)

list_to_std_vector_int()  
#list_to_std_list_int()   # it fails!
list_to_std_vector_dd()   # it fails! 
