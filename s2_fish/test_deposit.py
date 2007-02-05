#!/usr/bin/env python

from s2_fish import *
from macro_bunch import Macro_bunch
import Numeric

import unittest
class Test_deposit(unittest.TestCase):
    def test_01_ngp1(self):
        num_grid = 5
        num_grid = 5
        sf = Real_scalar_field(int3(num_grid,num_grid,num_grid),
                          double3(2.0,2.0,2.0),double3(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=1.0)
        total_charge = deposit_charge_ngp(sf,mb.store)
        self.assertAlmostEqual(total_charge,mb.store.local_num)
        self.assertAlmostEqual(sf.get_point(int3(0,2,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(4,2,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,0,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,4,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,2,0)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,2,4)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,2,2)),512.0)
        self.assertAlmostEqual(sf.get_point(int3(3,3,1)),64.0)
        self.assertAlmostEqual(sf.get_point(int3(2,1,3)),128.0)

    def test_02_ngp2(self):
        num_grid = 64
        sf = Real_scalar_field(int3(num_grid,num_grid,num_grid),
                          double3(2.0,2.0,2.0),double3(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=4.0)
        total_charge = deposit_charge_ngp(sf,mb.store)
        self.assertAlmostEqual(total_charge/mb.store.local_num,1.0/8.0)
        
    def test_03_cic1(self):
        num_grid = 5
        sf = Real_scalar_field(int3(num_grid,num_grid,num_grid),
                          double3(2.0,2.0,2.0),double3(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=1.0)
        total_charge = deposit_charge_cic(sf,mb.store)
        self.assertAlmostEqual(total_charge,mb.store.local_num)
        self.assertAlmostEqual(sf.get_point(int3(0,2,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(4,2,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,0,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,4,2)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,2,0)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,2,4)),0.0)
        self.assertAlmostEqual(sf.get_point(int3(2,2,2)),512.0)
        self.assertAlmostEqual(sf.get_point(int3(3,3,1)),64.0)
        self.assertAlmostEqual(sf.get_point(int3(2,1,3)),128.0)

    def test_04_cic2(self):
        num_grid = 64
        sf = Real_scalar_field(int3(num_grid,num_grid,num_grid),
                          double3(2.0,2.0,2.0),double3(0.0,0.0,0.0))
        mb = Macro_bunch()
        mb.init_test(16,edge_length=4.0)
        total_charge = deposit_charge_cic(sf,mb.store)
        self.assertAlmostEqual(total_charge/mb.store.local_num,1.0/8.0)

if __name__ == '__main__':
    deposit_suite = unittest.TestLoader().loadTestsFromTestCase(Test_deposit)
    unittest.TextTestRunner(verbosity=2).run(deposit_suite)
