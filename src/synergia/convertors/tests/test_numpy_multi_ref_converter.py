#!/usr/bin/env python

import sys
#sys.path.append('..')

import test_helper
import numpy
from nose.tools import *

def f(i, j=0, k=0):
    return i + 10*j + 100*k + 0.25

n1 = 2
n2 = 3
n3 = 4

def test_convert_MArray1d():
    ah = test_helper.Array_holder(n1, n2, n3)
    for i in range(0, n1):
        assert_equal(ah.get_1d()[i], f(i))

def test_convert_MArray2d():
    ah = test_helper.Array_holder(n1, n2, n3)
    for i in range(0, n1):
        for j in range(0, n2):
            assert_equal(ah.get_2d()[i][j], f(i, j))

def test_convert_MArray3d():
    ah = test_helper.Array_holder(n1, n2, n3)
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                assert_equal(ah.get_3d()[i][j][k], f(i, j, k))

def test_convert_MArray2d_fortran():
    ah = test_helper.Array_holder(n1, n2, n3)
    for i in range(0, n1):
        for j in range(0, n2):
            assert_equal(ah.get_2d_fortran()[i][j], f(i, j))

def test_convert_MArray3d_fortran():
    ah = test_helper.Array_holder(n1, n2, n3)
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                assert_equal(ah.get_3d_fortran()[i][j][k], f(i, j, k))

