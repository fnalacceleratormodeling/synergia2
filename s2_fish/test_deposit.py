#!/usr/bin/env python

from s2_fish import *
from macro_bunch import Macro_bunch
import Numeric
import sys

import unittest
class Test_deposit(unittest.TestCase):
    def test_01_ngp1(self):
        num_grid = 5
        sf = Real_scalar_field((num_grid,num_grid,num_grid),
                          (2.0,2.0,2.0),(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=1.0)
        h = sf.get_cell_size()
        vol = h[0]*h[1]*h[2]
        total_charge = deposit_charge_ngp(sf,mb.get_store())
        self.assertAlmostEqual(total_charge,mb.get_store().local_num)
        self.assertAlmostEqual(sf.get_points().get((0,2,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((4,2,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,0,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,4,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,2,0)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,2,4)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,2,2)),512.0/vol)
        self.assertAlmostEqual(sf.get_points().get((3,3,1)),64.0/vol)
        self.assertAlmostEqual(sf.get_points().get((2,1,3)),128.0/vol)

    def test_02_ngp2(self):
        num_grid = 64
        sf = Real_scalar_field((num_grid,num_grid,num_grid),
                          (2.0,2.0,2.0),(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=4.0)
        total_charge = deposit_charge_ngp(sf,mb.get_store())
        self.assertAlmostEqual(total_charge/mb.get_store().local_num,1.0/8.0)
        
    def test_03_cic1(self):
        num_grid = 5
        sf = Real_scalar_field((num_grid,num_grid,num_grid),
                          (2.0,2.0,2.0),(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=1.0)
        h = sf.get_cell_size()
        vol = h[0]*h[1]*h[2]
        total_charge = deposit_charge_cic(sf,mb.get_store(),0)
        self.assertAlmostEqual(total_charge,mb.get_store().local_num)
        self.assertAlmostEqual(sf.get_points().get((0,2,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((4,2,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,0,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,4,2)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,2,0)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,2,4)),0.0)
        self.assertAlmostEqual(sf.get_points().get((2,2,2)),512.0/vol)
        self.assertAlmostEqual(sf.get_points().get((3,3,1)),64.0/vol)
        self.assertAlmostEqual(sf.get_points().get((2,1,3)),128.0/vol)

    def test_04_cic2(self):
        num_grid = 128
        physical_size = (10.0,10.0,10.0)
        sf = Real_scalar_field((num_grid,num_grid,num_grid),
                          physical_size,(0.0,0.0,0.0))
        mb = Macro_bunch()
        num_per_side = 10
        mb.init_test(num_per_side,edge_length=4.0)
        h = sf.get_cell_size()
        vol = h[0]*h[1]*h[2]
        total_charge = deposit_charge_cic(sf,mb.get_store(),0)
        self.assertAlmostEqual(total_charge,num_per_side**3)

if __name__ == '__main__':
    unsuccessful = 0
    deposit_suite = unittest.TestLoader().loadTestsFromTestCase(Test_deposit)
    retval = unittest.TextTestRunner(verbosity=2).run(deposit_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
