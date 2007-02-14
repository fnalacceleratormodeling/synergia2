#!/usr/bin/env python

import sys

from s2_fish import *
import unittest
class Test_Real_scalar_field(unittest.TestCase):
    def test_01_construct1(self):
        sf = Real_scalar_field()

    def test_02_contruct2(self):
        sf = Real_scalar_field((10,10,10),(1.0,1.0,1.0),
                          (0.0,0.0,0.0))

    def test_03_num_points(self):
        sf = Real_scalar_field((7,17,37),(1.0,1.0,1.0),
                          (0.0,0.0,0.0))
        num = sf.get_points().get_shape()
        self.assertEqual(num[0],7)
        self.assertEqual(num[1],17)
        self.assertEqual(num[2],37)

    def test_04_physical_size(self):
        sf = Real_scalar_field((10,10,10),(1.0,2.0,3.0),
                          (0.0,0.0,0.0))
        size = sf.get_physical_size()
        self.assertEqual(size[0],1.0)
        self.assertEqual(size[1],2.0)
        self.assertEqual(size[2],3.0)

    def test_05_physical_offset(self):
        sf = Real_scalar_field((10,10,10),(1.0,2.0,3.0),
                          (0.1,0.2,0.3))
        offset = sf.get_physical_offset()
        self.assertEqual(offset[0],0.1)
        self.assertEqual(offset[1],0.2)
        self.assertEqual(offset[2],0.3)

    def test_06_set_get(self):
        sf = Real_scalar_field((10,10,10),(1.0,2.0,3.0),
                          (0.1,0.2,0.3))
        sqrtpi = 1.7724539
        sf.get_points().set((3,5,7),sqrtpi)
        val = sf.get_points().get((3,5,7))
        self.assertEqual(val,sqrtpi)
        
    def test_07_set_add_get(self):
        sf = Real_scalar_field((10,10,10),(1.0,2.0,3.0),
                          (0.1,0.2,0.3))
        sqrtpi = 1.7724539
        e = 2.718281828459045
        sf.get_points().set((4,6,8),sqrtpi)
        sf.get_points().add_to_point((4,6,8),e)
        val = sf.get_points().get((4,6,8))
        self.assertEqual(val,sqrtpi+e)

    def test_08_set_zero_get(self):
        sf = Real_scalar_field((10,10,10),(1.0,2.0,3.0),
                          (0.1,0.2,0.3))
        sqrtpi = 1.7724539
        sf.get_points().set((4,5,3),sqrtpi)
        sf.get_points().zero_all()
        val = sf.get_points().get((4,5,3))
        self.assertEqual(val,0.0)

    def test_09_set_outside(self):
        sf = Real_scalar_field((10,10,10),(1.0,2.0,3.0),
                          (0.1,0.2,0.3))
        caught = 0
        try:
            sf.get_points().set((10,9,9),1.0)
        except IndexError,e:
            caught = 1
        self.assertEqual(caught,1)

    def test_10_get_outside(self):
        sf = Real_scalar_field((10,10,10),(1.0,2.0,3.0),
                          (0.1,0.2,0.3))
        caught = 0
        try:
            retval = sf.get_points().get((10,9,9))
        except IndexError,e:
            caught = 1
        self.assertEqual(caught,1)

    def test_11_get_leftmost_indices(self):
        n = (10,10,10)
        size = (1.0,2.0,3.0)
        offset = (0.1,0.2,0.3)
        testpoint = (0.44,0.99,-1.01)
        h = (size[0]/(n[0]-1),
             size[1]/(n[1]-1),
             size[2]/(n[2]-1))
        left = (offset[0] - 0.5*size[0],
                offset[1] - 0.5*size[1],
                offset[2] - 0.5*size[2])
        location = ((testpoint[0] - left[0])/h[0],
                    (testpoint[1] - left[1])/h[1],
                    (testpoint[2] - left[2])/h[2])
        expected = tuple(map(int,location))
        
        sf = Real_scalar_field(n,size,offset)
        retval = sf.get_leftmost_indices(testpoint)
        self.assertEqual(retval,expected)

    def test_12_get_leftmost_offsets(self):
        sf = Real_scalar_field((2,2,2),(10.0,10.0,10.0),
                          (5.0,5.0,5.0))
        retval = sf.get_leftmost_offsets((1.0,3.0,5.0))
        self.assertAlmostEqual(retval[0],0.1,12)
        self.assertAlmostEqual(retval[1],0.3,12)
        self.assertAlmostEqual(retval[2],0.5,12)
        

if __name__ == '__main__':
    unsuccessful = 0
    scalar_field_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Real_scalar_field)
    retval = unittest.TextTestRunner(verbosity=2).run(scalar_field_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
