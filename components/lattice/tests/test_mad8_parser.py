#!/usr/bin/env python

import sys
sys.path.append('..')

from nose.tools import *
from mad8_parser import Expression_parser
    
    
#    ep = Expression_parser()
##    stack = ep.parse('3*4+5*pi/2.0')
#    stack = ep.parse('2*pi+4*5')
#    print 'ans =', ep.evaluate_stack(stack)
#
#    stack = ep.parse('sin(2.3)')
#    print 'ans =', ep.evaluate_stack(stack)
#
##    stack = ep.parse('-sin(2.3)')
##    print 'ans =', ep.evaluate_stack(stack)
#
#    stack = ep.parse('cos(x)')
#    print 'ans =', ep.evaluate_stack(stack)
    
def test_construct():
    ep = Expression_parser()

def test_integer():
    ep = Expression_parser()
    num = 1234
    stack = ep.parse(str(num))
    assert_equal(num,ep.evaluate_stack(stack))

def test_neg_integer():
    ep = Expression_parser()
    num = -3456
    stack = ep.parse(str(num))
    assert_equal(num,ep.evaluate_stack(stack))
    
def test_float():
    ep = Expression_parser()
    num = 12.34
    stack = ep.parse(str(num))
    assert_equal(num,ep.evaluate_stack(stack))

def test_neg_float():
    ep = Expression_parser()
    num = -12.34
    stack = ep.parse(str(num))
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_leading_dot():
    ep = Expression_parser()
    num = 0.678
    stack = ep.parse('.678')
    assert_equal(num,ep.evaluate_stack(stack))

def test_neg_float_leading_dot():
    ep = Expression_parser()
    num = -0.678
    stack = ep.parse('-.678')
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_e():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35e19')
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_E():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35E19')
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_D():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35D19')
    assert_equal(num,ep.evaluate_stack(stack))

