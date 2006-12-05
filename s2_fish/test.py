#!/usr/bin/env python

from s2_fish import *
import unittest
class Test_Scalar_Field(unittest.TestCase):
    def test_01_construct1(self):
        sf = Scalar_Field()

    def test_02_contruct2(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,1.0,1.0),
                          double3(0.0,0.0,0.0))

    def test_03_num_points1(self):
        sf = Scalar_Field()
        sf.set_num_points(int3(7,17,37))
        num = sf.get_num_points()
        self.assertEqual(num.get(0),7)
        self.assertEqual(num.get(1),17)
        self.assertEqual(num.get(2),37)

if __name__ == '__main__':
    scalar_field_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Scalar_Field)
    unittest.TextTestRunner(verbosity=2).run(scalar_field_suite)
