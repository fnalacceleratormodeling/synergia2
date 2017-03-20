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
val = 3.14

def test_convert_MArray1d():
    a = test_helper.get_MArray1d(n1)
    for i in range(0, n1):
        assert_equal(a[i], f(i))
    a[0] = val
    assert_equal(a[0], val)


def test_convert_MArray2d():
    a = test_helper.get_MArray2d(n1, n2)
    for i in range(0, n1):
        for j in range(0, n2):
            assert_equal(a[i][j], f(i, j))
    a[0][0] = val
    assert_equal(a[0][0], val)

def test_convert_MArray3d():
    a = test_helper.get_MArray3d(n1, n2, n3)
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                assert_equal(a[i][j][k], f(i, j, k))
    a[0][0][0] = val
    assert_equal(a[0][0][0], val)

def test_convert_MArray2d_fortran():
    a = test_helper.get_MArray2d_fortran(n1, n2)
    for i in range(0, n1):
        for j in range(0, n2):
            assert_equal(a[i][j], f(i, j))
    a[0][0] = val
    assert_equal(a[0][0], val)

def test_convert_MArray3d_fortran():
    a = test_helper.get_MArray3d_fortran(n1, n2, n3)
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                assert_equal(a[i][j][k], f(i, j, k))
    a[0][0][0] = val
    assert_equal(a[0][0][0], val)
