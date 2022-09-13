#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.append('..')
sys.path.append('../..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')

from mpi4py import MPI
from foundation import Reference_particle, Four_momentum
from bunch import Bunch, Bunch_train
from utils import Commxx, generate_subcomms
import convertors
import numpy
from nose.tools import *

mass = 2.2;
total_energy = 3.0
total_num = 100
real_num = 1.0e12
proton_charge = 1
num_bunches = 4
separation = 2.3
reference_particle = Reference_particle(proton_charge, mass, total_energy)

def test_construct1():
    comm=Commxx()
    commxxs = generate_subcomms(comm, num_bunches)
    bunches = []
    for commxx in commxxs:
        bunches.append(Bunch(reference_particle, total_num, real_num, commxx))
    bunch_train = Bunch_train(bunches, separation)

def test_construct2():
    comm=Commxx()
    commxxs = generate_subcomms(comm,num_bunches)
    bunches = []
    for commxx in commxxs:
        bunches.append(Bunch(reference_particle, total_num, real_num, commxx))
    separations = [separation, separation, separation]
    bunch_train = Bunch_train(bunches, separations)
