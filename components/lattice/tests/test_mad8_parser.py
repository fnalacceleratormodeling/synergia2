#!/usr/bin/env python

import sys
sys.path.append('..')

from nose.tools import *
from mad8_parser import Expression_parser
from math import sqrt,log,exp,sin,cos,tan,asin    
    
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
    
def test_pos_integer():
    ep = Expression_parser()
    num = 3456
    stack = ep.parse('+3456')
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

def test_float_e2():
    ep = Expression_parser()
    num = 1.35e-19
    stack = ep.parse('1.35e-19')
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_e3():
    ep = Expression_parser()
    num = 1.35e+19
    stack = ep.parse('1.35e+19')
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_E():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35E19')
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_d():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35d19')
    assert_equal(num,ep.evaluate_stack(stack))

def test_float_D():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35D19')
    assert_equal(num,ep.evaluate_stack(stack))

def test_simple_expression():
    ep = Expression_parser()
    stack = ep.parse('1.1+2.2*(3.3+4.4)^6.6')
    assert_almost_equal(1.1+2.2*(3.3+4.4)**6.6,ep.evaluate_stack(stack))

def test_functions():
    ep = Expression_parser()
    stack = ep.parse('sqrt(1.1)+log(2.2)+exp(3.3)+sin(4.4)+cos(5.5)')
    assert_almost_equal(sqrt(1.1)+log(2.2)+exp(3.3)+sin(4.4)+cos(5.5),
                        ep.evaluate_stack(stack))
    
def test_functions2():
    ep = Expression_parser()
    stack = ep.parse('tan(6.6)+asin(0.77)+abs(-8.8)')
    assert_almost_equal(tan(6.6)+asin(0.77)+abs(-8.8),
                        ep.evaluate_stack(stack))
    
    
