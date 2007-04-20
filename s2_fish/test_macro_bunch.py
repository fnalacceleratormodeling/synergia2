#!/usr/bin/env python

from macro_bunch import *

import sublocal_paths
import physics_constants
import beam_parameters
import bunch
import processor_grid

import math
import sys
import unittest

import syn2_diagnostics

class Test_Macro_bunch(unittest.TestCase):
    def test_01_construct(self):
        mb = Macro_bunch()

    def test_02_init_test1(self):
        mb = Macro_bunch()
        num_per_side = 2
        mb.init_test(num_per_side)
        shape = Numeric.shape(mb.get_local_particles())
        self.assertEqual(shape[0],7)
        self.assertEqual(shape[1],num_per_side**3)

    def test_03_init_test2(self):
        mb = Macro_bunch()
        num_per_side = 2
        mb.init_test(num_per_side)
        shape = Numeric.shape(mb.particles)
        for j in range(0,shape[1]):
            for i in range(0,shape[0]):
                self.assertAlmostEqual(mb.particles[i,j],
                                       mb.get_store().get_coord(i,j),14)
                
    def _get_bunch(self):
        current = 0.5
        kinetic_energy = 0.0067
        mass = physics_constants.PH_NORM_mp
        charge = 1.0
        initial_phase = 0.0
        scaling_frequency = 10221.05558e6
        width_x = 0.004
        num_particles = 1000
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

        return b
    
    def test_04_init_from_bunch(self):
        b = self._get_bunch()
        orig_shape = Numeric.shape(b.particles())
        mb = Macro_bunch()
        mb.init_from_bunch(b)
        shape = Numeric.shape(mb.get_local_particles())
        for j in range(0,shape[1]):
            for i in range(0,shape[0]):
                self.assertAlmostEqual(b.particles()[i,j],
                                       mb.get_store().get_coord(i,j),14)
        self.assertEqual(shape[0],orig_shape[0])
        self.assertEqual(shape[1],orig_shape[1])
        

    def test_05_conversions(self):
        b = self._get_bunch()
        mb = Macro_bunch()
        mb.init_from_bunch(b)
        sample_id = 307
        x0 = mb.get_local_particles()[0,sample_id]
        y0 = mb.get_local_particles()[2,sample_id]
        z0 = mb.get_local_particles()[4,sample_id]
        x1 = x0/mb.get_store().get_units()[0]
        y1 = y0/mb.get_store().get_units()[2]
        mb.convert_to_fixedt()
        self.assertEqual(x1,mb.get_local_particles()[0,sample_id])
        self.assertEqual(y1,mb.get_local_particles()[2,sample_id])
        # really need to check z conversion
        self.assertEqual(0,mb.get_store().is_fixedz)

        mb.convert_to_fixedz()
        self.assertEqual(x0,mb.get_local_particles()[0,sample_id])
        self.assertEqual(y0,mb.get_local_particles()[2,sample_id])
        self.assertEqual(z0,mb.get_local_particles()[4,sample_id])
        self.assertEqual(1,mb.get_store().is_fixedz)

    def test_06_diagnostics(self):
        b = self._get_bunch()
        orig_shape = Numeric.shape(b.particles())
        mb = Macro_bunch()
        mb.init_from_bunch(b)
        d = syn2_diagnostics.Diagnostics(mb.units)
        d.get_coord_stds(mb)

if __name__ == '__main__':
    unsuccessful = 0
    macro_bunch_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Macro_bunch)
    retval = unittest.TextTestRunner(verbosity=2).run(macro_bunch_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
