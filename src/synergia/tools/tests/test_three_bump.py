#!/usr/bin/env python

import sys
import os
import synergia
import synergia.tools
import numpy as np
import tables
import pygsl.errno as pygslerrno
import pygsl
from pygsl import multiroots
from nose.tools import *

#  just a little tester for the class
def test_three_bump():
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "lattices/foborodobo128.madx")
    print "read lattice: ", len(lattice.get_elements()), " elements, length = ", lattice.get_length()
    hcorr_names = ('hc1', 'hc2', 'hc3')
    vcorr_names = ('vc1', 'vc2', 'vc3')
    three_bump = synergia.tools.Three_bump(lattice, 'm1', 'm2', hcorr_names, vcorr_names, 'm3', False)
    three_bump.information()

    target_x = 0.001
    target_y = -0.0005
    bump_settings = three_bump.set_bump((target_x, target_y))
    print "bump_settings: ", bump_settings[0], bump_settings[1], bump_settings[2], bump_settings[3], bump_settings[4], bump_settings[5]

    # propagate the whole lattice now
    comm = synergia.utils.Commxx()
    refpart = lattice.get_reference_particle()
    stepper = synergia.simulation.Independent_stepper_elements(lattice, 1, 1)
    # 3 particles is the minimum so that the diagnostics don't crash
    bunch = synergia.bunch.Bunch(refpart, 3, 1.0e10, comm)
    bunch.get_local_particles()[:,0:6] = 0.0
    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
    bunch_simulator.add_per_step(synergia.bunch.Diagnostics_basic("step_basic.h5"))
    propagator = synergia.simulation.Propagator(stepper)
    propagator.propagate(bunch_simulator, 1, 1, 1)
    h5 = tables.openFile("bump_basic.h5")
    mean = h5.root.mean.read()
    h5.close()
    assert_almost_equal(mean[0,0], target_x)
    assert_almost_equal(mean[2,0], target_y)
    for i in range(4):
        assert_almost_equal(mean[i, 1], 0.0)
