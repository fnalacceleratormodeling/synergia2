#!/usr/bin/env python

from s2_fish import *
from macro_bunch import Macro_bunch
import Numeric
import time
import sys

import unittest
class Test_solver_fft_open(unittest.TestCase):    
    def test_01(self):
        max_err = fft_tester(8,32,64);
        self.assertAlmostEqual(max_err,0.0)
        
    def test_02(self):
        num_grid = 32
        sf = Real_scalar_field((32,32,32),(2.0,3.0,4.0),(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_sphere(100000,0.5)
        t0 = time.time()
        total_charge = deposit_charge_cic(sf,mb.store)
        t0 = time.time()
        phi = solver_fft_open(sf)

#     def test_03(self):
#         sf = Real_scalar_field()
#         sf.read_from_file("rho_in")
#         phi = solver(sf)

if __name__ == '__main__':
    unsuccessful = 0
    solver_suite = unittest.TestLoader().loadTestsFromTestCase(Test_solver_fft_open)
    retval = unittest.TextTestRunner(verbosity=2).run(solver_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
