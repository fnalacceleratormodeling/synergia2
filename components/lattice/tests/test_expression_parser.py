#!/usr/bin/env python

import sys
sys.path.append('..')

from nose.tools import *
from mad8_parser import Expression_parser, ParseException
from math import sqrt, log, exp, sin, cos, tan, asin    
import math

def test_construct():
    ep = Expression_parser()

def test_integer():
    ep = Expression_parser()
    num = 1234
    stack = ep.parse(str(num))
    assert_equal(num, ep.evaluate_stack(stack))

def test_neg_integer():
    ep = Expression_parser()
    num = - 3456
    stack = ep.parse(str(num))
    assert_equal(num, ep.evaluate_stack(stack))
    
def test_pos_integer():
    ep = Expression_parser()
    num = 3456
    stack = ep.parse('+3456')
    assert_equal(num, ep.evaluate_stack(stack))
    
def test_float():
    ep = Expression_parser()
    num = 12.34
    stack = ep.parse(str(num))
    assert_equal(num, ep.evaluate_stack(stack))

def test_neg_float():
    ep = Expression_parser()
    num = - 12.34
    stack = ep.parse(str(num))
    assert_equal(num, ep.evaluate_stack(stack))

def test_float_leading_dot():
    ep = Expression_parser()
    num = 0.678
    stack = ep.parse('.678')
    assert_equal(num, ep.evaluate_stack(stack))

def test_neg_float_leading_dot():
    ep = Expression_parser()
    num = - 0.678
    stack = ep.parse('-.678')
    assert_equal(num, ep.evaluate_stack(stack))

def test_float_e():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35e19')
    assert_equal(num, ep.evaluate_stack(stack))

def test_float_e2():
    ep = Expression_parser()
    num = 1.35e-19
    stack = ep.parse('1.35e-19')
    assert_equal(num, ep.evaluate_stack(stack))

def test_float_e3():
    ep = Expression_parser()
    num = 1.35e+19
    stack = ep.parse('1.35e+19')
    assert_equal(num, ep.evaluate_stack(stack))

def test_float_E():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35E19')
    assert_equal(num, ep.evaluate_stack(stack))

def test_float_d():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35d19')
    assert_equal(num, ep.evaluate_stack(stack))

def test_float_D():
    ep = Expression_parser()
    num = 1.35e19
    stack = ep.parse('1.35D19')
    assert_equal(num, ep.evaluate_stack(stack))

def test_simple_expression():
    ep = Expression_parser()
    stack = ep.parse('1.1+2.2*(3.3+4.4)^6.6')
    assert_almost_equal(1.1 + 2.2 * (3.3 + 4.4) ** 6.6, ep.evaluate_stack(stack))

def test_functions():
    ep = Expression_parser()
    stack = ep.parse('sqrt(1.1)+log(2.2)+exp(3.3)+sin(4.4)+cos(5.5)')
    assert_almost_equal(sqrt(1.1) + log(2.2) + exp(3.3) + sin(4.4) + cos(5.5),
                        ep.evaluate_stack(stack))
    
def test_functions2():
    ep = Expression_parser()
    stack = ep.parse('tan(6.6)+asin(0.77)+abs(-8.8)')
    assert_almost_equal(tan(6.6) + asin(0.77) + abs(-8.8),
                        ep.evaluate_stack(stack))

def test_max():
    ep = Expression_parser()
    stack = ep.parse('max(-4e3,7)')
    assert_equal(7, ep.evaluate_stack(stack))
    
def test_min():
    ep = Expression_parser()
    stack = ep.parse('min(3,10.8)')
    assert_equal(3, ep.evaluate_stack(stack))
    
def test_function_neg():
    ep = Expression_parser()
    stack = ep.parse('-sin(2.3)')
    assert_almost_equal(-sin(2.3), ep.evaluate_stack(stack))
   
def test_pi():
    ep = Expression_parser()
    stack = ep.parse('pi')
    assert_almost_equal(math.pi, ep.evaluate_stack(stack))

def test_twopi():
    ep = Expression_parser()
    stack = ep.parse('twopi')
    assert_almost_equal(2 * math.pi, ep.evaluate_stack(stack))

def test_degrad():
    ep = Expression_parser()
    stack = ep.parse('degrad')
    assert_almost_equal(180.0 / math.pi, ep.evaluate_stack(stack))

def test_raddeg():
    ep = Expression_parser()
    stack = ep.parse('raddeg')
    assert_almost_equal(math.pi / 180.0, ep.evaluate_stack(stack))

def test_e():
    ep = Expression_parser()
    stack = ep.parse('e')
    assert_almost_equal(math.e, ep.evaluate_stack(stack))

def test_emass():
    ep = Expression_parser()
    stack = ep.parse('emass')
    assert_almost_equal(0.51099906e-3, ep.evaluate_stack(stack))

def test_pmass():
    ep = Expression_parser()
    stack = ep.parse('pmass')
    assert_almost_equal(0.93827231, ep.evaluate_stack(stack))

def test_clight():
    ep = Expression_parser()
    stack = ep.parse('clight')
    assert_almost_equal(2.99792458e8, ep.evaluate_stack(stack))

def test_vars():
    ep = Expression_parser()
    a = 1.1
    b = 2.2
    vars = {'a':a, 'b':b}
    stack = ep.parse('a+b')
    assert_almost_equal(a + b, ep.evaluate_stack(stack, vars))

def test_vars2():
    ep = Expression_parser()
    a = 1.1
    b = 2.2
    c = 3.3
    d = 4.4
    vars = {'a':a, 'b':b, 'c':c, 'd':d}
    stack = ep.parse('a+b*(c+d)')
    assert_almost_equal(a + b * (c + d), ep.evaluate_stack(stack, vars))

def test_unassigned_var():
    ep = Expression_parser(uninitialized_warn=False)
    stack = ep.parse('foo')
    assert_equal(0.0, ep.evaluate_stack(stack))

def test_unassigned_var_exception():
    ep = Expression_parser(uninitialized_exception=True)
    stack = ep.parse('foo')
    caught = False
    try:
        value = ep.evaluate_stack(stack)
    except ParseException, e:
        caught = True
    assert caught

