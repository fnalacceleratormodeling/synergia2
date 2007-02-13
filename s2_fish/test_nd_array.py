#!/usr/bin/env python

from nd_array import Real_nd_array, Complex_nd_array
import tempfile,os
import unittest
class Test_Real_nd_array(unittest.TestCase):
    def create222(self):
        na = Real_nd_array([2,2,2])
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    na.set([i,j,k],i*100+j*10+k)
        return na

    def test_01_construct1(self):
        na = Real_nd_array()

    def test_02_construct2(self):
        na = Real_nd_array([2,2])

    def test_03_get_shape(self):
        shape_in = (2,3)
        na = Real_nd_array(shape_in)
        shape_out = na.get_shape()
        self.assertEqual(shape_in,shape_out)

    def test_04_copy(self):
        na_orig = self.create222()
        na_new = Real_nd_array()
        na_new.copy(na_orig)
        self.assertEqual(na_orig.get_shape(),na_new.get_shape())
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    self.assertEqual(na_orig.get([i,j,k]),
                                     na_new.get([i,j,k]))
                    
    def test_05_reshape(self):
        na = Real_nd_array()
        shape_out = na.get_shape()
        self.assertEqual(len(shape_out),0)
        shape_in = (2,3)
        na.reshape(shape_in) 
        shape_out = na.get_shape()
        self.assertEqual(shape_in,shape_out)

    def test_06_set_nocheck(self):
        na = Real_nd_array([5,5])
        na.set_nocheck([1,3],100.0)
        self.assertAlmostEqual(na.get([1,3]),100.0)        

    def test_07_zero_all(self):
        na = Real_nd_array([2,2])
        na.set_nocheck([0,0],100)
        na.zero_all()
        self.assertEqual(na.get([0,0]),0.0)

    def test_08_set(self):
        na = Real_nd_array([5,5,7])
        na.set([1,3,6],100.0)
        self.assertAlmostEqual(na.get([1,3,6]),100.0)        

    def test_09_set_out_of_bounds1(self):
        na = Real_nd_array([3,5,7,9])
        caught_index = 0
        try:
            na.set([5,0,0,0],100.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_10_set_out_of_bounds2(self):
        na = Real_nd_array([3,5,7,9])
        caught_index = 0
        try:
            na.set([2,4,6,9],100.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_11_set_out_of_bounds3(self):
        na = Real_nd_array([3,5,7,9])
        caught_index = 0
        try:
            na.set([0,0,-1,0],100.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_12_set_out_of_bounds4(self):
        na = Real_nd_array([16,32,8])
        caught_num = 0
        try:
            na.set([16,0,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,32,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,0,8],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([-1,0,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,-1,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,0,-1],100.0)
        except IndexError,e:
            caught_num += 1
        self.assertEqual(caught_num,6)

    def test_13_set_empty(self):
        na = Real_nd_array()
        caught_index = 0
        try:
            na.set([0],0.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_14_scale(self):
        na = Real_nd_array([2,2])
        na.set([0,0],1.0)
        na.set([0,1],2.0)
        na.set([1,0],3.0)
        na.set([1,1],4.0)
        na.scale(101.0)
        self.assertEqual(na.get([0,0]),101.0)
        self.assertEqual(na.get([0,1]),202.0)
        self.assertEqual(na.get([1,0]),303.0)
        self.assertEqual(na.get([1,1]),404.0)

    def test_15_get_length(self):
        shape = (3,5,17)
        expected_length = shape[0]*shape[1]*shape[2]
        na = Real_nd_array(shape)
        length = na.get_length()
        self.assertEqual(length,expected_length)

    def test_16_get_length2(self):
        na = Real_nd_array()        
        length = na.get_length()
        self.assertEqual(length,0)
        shape = (2,64,8,5)
        expected_length = shape[0]*shape[1]*shape[2]*shape[3]
        na.reshape(shape)
        length = na.get_length()
        self.assertEqual(length,expected_length)

    def test_17_describe(self):
        print
        na = Real_nd_array()
        na.describe()
        na = Real_nd_array([1,2,3,4,5,6,7])
        na.describe()

    def test_18_print_empty(self):
        print
        empty = Real_nd_array()
        empty.print_("empty")

    def test_19_print(self):
        print
        na = self.create222()
        na.print_("na")

    def test_20_print(self):
        print
        na = Real_nd_array([4])
        na.print_("na")

    def test_21_readwrite(self):
        na_orig = self.create222()
        filename = tempfile.mktemp()
        na_orig.write_to_file(filename)
        na_new = Real_nd_array()
        na_new.read_from_file(filename)
        self.assertEqual(na_orig.get_shape(),na_new.get_shape())
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    self.assertEqual(na_orig.get([i,j,k]),
                                     na_new.get([i,j,k]))
        os.unlink(filename)

class Test_Complex_nd_array(unittest.TestCase):
    def create222(self):
        na = Complex_nd_array([2,2,2])
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    na.set([i,j,k],i*10.0+j*1.0+k*1.0j)
        return na

    def test_01_construct1(self):
        na = Complex_nd_array()

    def test_02_construct2(self):
        na = Complex_nd_array([2,2])

    def test_03_get_shape(self):
        shape_in = (2,3)
        na = Complex_nd_array(shape_in)
        shape_out = na.get_shape()
        self.assertEqual(shape_in,shape_out)

    def test_04_copy(self):
        na_orig = self.create222()
        na_new = Complex_nd_array()
        na_new.copy(na_orig)
        self.assertEqual(na_orig.get_shape(),na_new.get_shape())
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    self.assertEqual(na_orig.get([i,j,k]),
                                     na_new.get([i,j,k]))

    def test_05_copy_real(self):
        na_real = Real_nd_array([2,2,2])
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    na_real.set([i,j,k],i+100+j+10+k)
        na = Complex_nd_array()
        na.copy_real(na_real)
        self.assertEqual(na_real.get_shape(),na.get_shape())
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    self.assertEqual(na_real.get([i,j,k]),
                                     na.get([i,j,k]).real)                    
                    self.assertEqual(0.0,
                                     na.get([i,j,k]).imag)
                    
    def test_06_reshape(self):
        na = Complex_nd_array()
        shape_out = na.get_shape()
        self.assertEqual(len(shape_out),0)
        shape_in = (2,3)
        na.reshape(shape_in) 
        shape_out = na.get_shape()
        self.assertEqual(shape_in,shape_out)

    def test_07_set_nocheck(self):
        na = Complex_nd_array([5,5])
        na.set_nocheck([1,3],100.0j)
        self.assertAlmostEqual(na.get([1,3]).real,0.0)        
        self.assertAlmostEqual(na.get([1,3]).imag,100.0)        

    def test_08_zero_all(self):
        na = Complex_nd_array([2,2])
        na.set_nocheck([0,0],100)
        na.zero_all()
        self.assertEqual(na.get([0,0]),0.0)

    def test_09_set(self):
        na = Complex_nd_array([5,5,7])
        na.set([1,3,6],100.0j)
        self.assertAlmostEqual(na.get([1,3,6]).real,0.0)        
        self.assertAlmostEqual(na.get([1,3,6]).imag,100.0)        

    def test_10_set_out_of_bounds1(self):
        na = Complex_nd_array([3,5,7,9])
        caught_index = 0
        try:
            na.set([5,0,0,0],100.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_11_set_out_of_bounds2(self):
        na = Complex_nd_array([3,5,7,9])
        caught_index = 0
        try:
            na.set([2,4,6,9],100.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_12_set_out_of_bounds3(self):
        na = Complex_nd_array([3,5,7,9])
        caught_index = 0
        try:
            na.set([0,0,-1,0],100.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_13_set_out_of_bounds4(self):
        na = Complex_nd_array([16,32,8])
        caught_num = 0
        try:
            na.set([16,0,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,32,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,0,8],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([-1,0,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,-1,0],100.0)
        except IndexError,e:
            caught_num += 1
        try:
            na.set([0,0,-1],100.0)
        except IndexError,e:
            caught_num += 1
        self.assertEqual(caught_num,6)

    def test_14_set_empty(self):
        na = Complex_nd_array()
        caught_index = 0
        try:
            na.set([0],0.0)
        except IndexError,e:
            caught_index = 1
        self.assertEqual(caught_index,1)

    def test_15_scale(self):
        na = Complex_nd_array([2,2])
        na.set([0,0],1.0)
        na.set([0,1],2.0)
        na.set([1,0],3.0)
        na.set([1,1],4.0)
        na.scale(101.0)
        self.assertEqual(na.get([0,0]),101.0)
        self.assertEqual(na.get([0,1]),202.0)
        self.assertEqual(na.get([1,0]),303.0)
        self.assertEqual(na.get([1,1]),404.0)

    def test_16_get_length(self):
        shape = (3,5,17)
        expected_length = shape[0]*shape[1]*shape[2]
        na = Complex_nd_array(shape)
        length = na.get_length()
        self.assertEqual(length,expected_length)

    def test_17_get_length2(self):
        na = Complex_nd_array()        
        length = na.get_length()
        self.assertEqual(length,0)
        shape = (2,64,8,5)
        expected_length = shape[0]*shape[1]*shape[2]*shape[3]
        na.reshape(shape)
        length = na.get_length()
        self.assertEqual(length,expected_length)

    def test_18_describe(self):
        print
        na = Complex_nd_array()
        na.describe()
        na = Complex_nd_array([1,2,3,4,5,6,7])
        na.describe()

    def test_19_print_empty(self):
        print
        empty = Complex_nd_array()
        empty.print_("empty")

    def test_20_print(self):
        print
        na = self.create222()
        na.print_("na")

    def test_21_print(self):
        print
        na = Complex_nd_array([4])
        na.print_("na")

    def test_22_readwrite(self):
        na_orig = self.create222()
        filename = tempfile.mktemp()
        na_orig.write_to_file(filename)
        na_new = Complex_nd_array()
        na_new.read_from_file(filename)
        self.assertEqual(na_orig.get_shape(),na_new.get_shape())
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    self.assertEqual(na_orig.get([i,j,k]),
                                     na_new.get([i,j,k]))
        os.unlink(filename)

if __name__ == '__main__':
    real_nd_array_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Real_nd_array)
    unittest.TextTestRunner(verbosity=2).run(real_nd_array_suite)
    complex_nd_array_suite = unittest.TestLoader().loadTestsFromTestCase(Test_Complex_nd_array)
    unittest.TextTestRunner(verbosity=2).run(complex_nd_array_suite)
