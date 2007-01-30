#!/usr/bin/env python

from s2_fish import *
from macro_bunch import Macro_bunch
import Numeric

import unittest
class Test_solver(unittest.TestCase):
    def test_01(self):
        num_grid = 8
        sf = Scalar_Field(int3(num_grid,num_grid,num_grid),
                          double3(2.0,2.0,2.0),double3(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=1.0)
        total_charge = deposit_charge_cic(sf,mb.store)
        phi = solver(sf)
        phi.print_points()

if __name__ == '__main__':
    solver_suite = unittest.TestLoader().loadTestsFromTestCase(Test_solver)
    unittest.TextTestRunner(verbosity=2).run(solver_suite)
