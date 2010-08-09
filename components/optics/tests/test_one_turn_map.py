#!/usr/bin/env python

import sys
sys.path.append("../../simulation")
sys.path.append("../../lattice")
sys.path.append("../../foundation")
sys.path.append("..")
sys.path.append("/home/amundson/work/synergia2-old_devel_1_0/install/lib")

from one_turn_map import linear_one_turn_map
from mad8_reader import Mad8_reader
from pysimulation import Collective_operator, Lattice_simulator, \
    Split_operator_stepper

def test_linear_one_turn_map():
    num_steps = 1
    map_order = 2
    lattice = Mad8_reader().get_lattice("fodo", "../../lattice/tests/fodo.lat")
    space_charge = Collective_operator("dummy")
    lattice_simulator = Lattice_simulator(lattice, map_order)
    stepper = Split_operator_stepper(lattice_simulator, space_charge, num_steps)
    map = linear_one_turn_map(lattice_simulator)

