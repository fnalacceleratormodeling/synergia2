#!/usr/bin/env python

from s2_fish import *
import unittest
class Test_Scalar_Field(unittest.TestCase):
    def test_01_construct1(self):
        sf = Scalar_Field()

    def test_02_contruct2(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,1.0,1.0),
                          double3(0.0,0.0,0.0))

    def test_03_num_points(self):
        sf = Scalar_Field()
        sf.set_num_points(int3(7,17,37))
        num = sf.get_num_points()
        self.assertEqual(num.get(0),7)
        self.assertEqual(num.get(1),17)
        self.assertEqual(num.get(2),37)

    def test_04_physical_size(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.0,0.0,0.0))
        size = sf.get_physical_size()
        self.assertEqual(size.get(0),1.0)
        self.assertEqual(size.get(1),2.0)
        self.assertEqual(size.get(2),3.0)

    def test_05_physical_offset(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.1,0.2,0.3))
        offset = sf.get_physical_offset()
        self.assertEqual(offset.get(0),0.1)
        self.assertEqual(offset.get(1),0.2)
        self.assertEqual(offset.get(2),0.3)

    def test_06_set_get(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.1,0.2,0.3))
        sqrtpi = 1.7724539
        sf.set_point(int3(3,5,7),sqrtpi)
        val = sf.get_point(int3(3,5,7))
        self.assertEqual(val,sqrtpi)
        
    def test_07_set_add_get(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.1,0.2,0.3))
        sqrtpi = 1.7724539
        e = 2.718281828459045
        sf.set_point(int3(4,6,8),sqrtpi)
        sf.add_to_point(int3(4,6,8),e)
        val = sf.get_point(int3(4,6,8))
        self.assertEqual(val,sqrtpi+e)

    def test_08_set_zero_get(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.1,0.2,0.3))
        sqrtpi = 1.7724539
        sf.set_point(int3(4,5,3),sqrtpi)
        sf.zero_the_points()
        val = sf.get_point(int3(4,5,3))
        self.assertEqual(val,0.0)

    def test_09_set_outside(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.1,0.2,0.3))
        caught = 0
        try:
            sf.set_point(int3(10,9,9),1.0)
        except IndexError,e:
            caught = 1
        self.assertEqual(caught,1)

    def test_10_get_outside(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.1,0.2,0.3))
        caught = 0
        try:
            retval = sf.get_point(int3(10,9,9))
        except IndexError,e:
            caught = 1
        self.assertEqual(caught,1)

    def test_11_get_leftmost_indices(self):
        sf = Scalar_Field(int3(10,10,10),double3(1.0,2.0,3.0),
                          double3(0.1,0.2,0.3))
        retval = sf.get_leftmost_indices(double3(0.44,0.99,-1.01))
        self.assertEqual(retval.get(0),7)
        self.assertEqual(retval.get(1),8)
        self.assertEqual(retval.get(2),0)

    def test_12_get_leftmost_offsets(self):
        sf = Scalar_Field(int3(2,2,2),double3(10.0,10.0,10.0),
                          double3(5.0,5.0,5.0))
        retval = sf.get_leftmost_offsets(double3(1.0,3.0,5.0))
        self.assertAlmostEqual(retval.get(0),0.1,12)
        self.assertAlmostEqual(retval.get(1),0.3,12)
        self.assertAlmostEqual(retval.get(2),0.5,12)
        

if __name__ == '__main__':
    scalar_field_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Scalar_Field)
    unittest.TextTestRunner(verbosity=2).run(scalar_field_suite)
