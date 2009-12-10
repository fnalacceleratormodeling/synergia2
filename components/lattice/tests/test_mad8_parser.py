#!/usr/bin/env python

import sys
sys.path.append('..')

from nose.tools import *
from mad8_parser import Mad8_parser, ParseException
from math import sqrt, log, exp, sin, cos, tan, asin    
import math

def test_construct():
    mp = Mad8_parser()

def test_blank_line():
    mp = Mad8_parser()
    mp.parse('')

def test_garbage():
    mp = Mad8_parser()
    caught = False
    try:
        mp.parse('your momma uses portions->of madx syntax')
    except ParseException, e:
        caught = True
    assert caught

def test_comment():    
    mp = Mad8_parser()
    mp.parse('!your momma uses portions->of madx syntax')

def test_variable_assignment():
    mp = Mad8_parser()
    mp.parse('x=1')
    assert_equal(1,mp.variables['x'])
    
def test_mod_variable_assignment():
    mp = Mad8_parser()
    mp.parse('x:=1')
    assert_equal(1,mp.variables['x'])
    
def test_newline_separation():
    mp = Mad8_parser()
    mp.parse('''x=1
    y=2''')
    assert_equal(1,mp.variables['x'])
    assert_equal(2,mp.variables['y'])
    
def test_semicolon_separation():
    mp = Mad8_parser()
    mp.parse('''x=1;y=2''')
    assert_equal(1,mp.variables['x'])
    assert_equal(2,mp.variables['y'])
    
def test_variable_assignment_expression():
    mp = Mad8_parser()
    mp.parse('foo.bar=pi*sin(1.2d-4)^0.69')
    assert_equal(math.pi*math.sin(1.2e-4)**0.69,mp.variables['foo.bar'])

def test_command():
    mp = Mad8_parser()
    mp.parse('foo')
    assert_equal(1,len(mp.commands))
    command = mp.commands[0]
    assert_equal('foo',command[0])
    assert_equal({},command[1])

def test_upper_command():
    mp = Mad8_parser()
    mp.parse('FOO')
    assert_equal(1,len(mp.commands))
    command = mp.commands[0]
    assert_equal('foo',command[0])
    assert_equal({},command[1])

def test_command_attrs():
    mp = Mad8_parser()
    mp.parse('foo,a=1,B=3*(4+5)')
    assert_equal(1,len(mp.commands))
    command = mp.commands[0]
    assert_equal('foo',command[0])
    attributes = command[1]
    assert_equal(2,len(attributes))
    assert_equal(1,attributes['a'])
    assert_equal(3*(4+5),attributes['b'])
    