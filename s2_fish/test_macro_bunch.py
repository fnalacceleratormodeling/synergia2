#!/usr/bin/env python

from macro_bunch import *

import sublocal_paths
import physics_constants
import beam_parameters
import bunch
import processor_grid

import math
import unittest

class Test_Macro_bunch(unittest.TestCase):
    def test_01_construct(self):
        mb = Macro_bunch()

    def test_02_init_test(self):
        mb = Macro_bunch()
        mb.init_test()

    def test_03_init_from_bunch(self):
        current = 0.5
        kinetic_energy = 0.0067
        mass = physics_constants.PH_NORM_mp
        charge = 1.0
        initial_phase = 0.0
        scaling_frequency = 10221.05558e6
        width_x = 0.004
        num_particles = 10000
        xwidth=0.0012026
        xpwidth=0.0049608
        rx=0.85440
        dpop = 1.0e-20

        bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                             initial_phase, scaling_frequency,
                                             transverse=1)

        pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
        bp.x_params(sigma = xwidth, lam = xpwidth * pz)
        bp.y_params(sigma = xwidth, lam = xpwidth * pz)
        sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c \
                         /scaling_frequency/math.pi
        bp.z_params(sigma = sigma_z_meters, lam = dpop* pz)
        bp.correlation_coeffs(xpx = -rx, ypy = rx)
        pgrid = processor_grid.Processor_grid(1)
        b = bunch.Bunch(current, bp, num_particles, pgrid)
        b.generate_particles()

        mb = Macro_bunch()
        mb.init_from_bunch(b)
        shape = Numeric.shape(mb.store.get_local_particles())
        self.assertEqual(shape[0],7)
        self.assertEqual(shape[1],10000)
        


if __name__ == '__main__':
    macro_bunch_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Macro_bunch)
    unittest.TextTestRunner(verbosity=2).run(macro_bunch_suite)
