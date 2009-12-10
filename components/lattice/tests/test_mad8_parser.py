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
    assert_equal('foo',command.name)
    assert_equal({},command.attributes)

def test_upper_command():
    mp = Mad8_parser()
    mp.parse('FOO')
    assert_equal(1,len(mp.commands))
    command = mp.commands[0]
    assert_equal('foo',command.name)
    assert_equal({},command.attributes)

def test_command_attrs():
    mp = Mad8_parser()
    mp.parse('foo,a=1,B=3*(4+5)')
    assert_equal(1,len(mp.commands))
    command = mp.commands[0]
    assert_equal('foo',command.name)
    attributes = command.attributes
    assert_equal(2,len(attributes))
    assert_equal(1,attributes['a'])
    assert_equal(3*(4+5),attributes['b'])

def test_command_str_attrs():
    mp = Mad8_parser()
    mp.parse('TITLE, S = "Tevatron Collider Run II Lattice"')
    assert_equal(1,len(mp.commands))
    command = mp.commands[0]
    assert_equal('title',command.name)
    attributes = command.attributes
    assert_equal(1,len(attributes))
    assert_equal("Tevatron Collider Run II Lattice",attributes['s'])

def test_command_particle_attr():
    mp = Mad8_parser()
    mp.parse('beam, particle=proton')
    assert_equal(1,len(mp.commands))
    command = mp.commands[0]
    assert_equal('beam',command.name)
    attributes = command.attributes
    assert_equal(1,len(attributes))
    assert_equal('proton',attributes['particle'])
    
def test_command_assign():
    mp = Mad8_parser()
    mp.parse('q1: quadrupole,l=3.14')
    assert_equal(1,len(mp.labels))
    key = 'q1'
    command = mp.labels[key]
    assert_equal('quadrupole',command.name)
    attributes = command.attributes
    assert_equal(1,len(attributes))
    assert_almost_equal(3.14,attributes['l'])

def test_line():
    mp = Mad8_parser()
    mp.parse('''f:quad,l=1,k1=0.1
    o:drift,l=2
    d:quad,l=1,k1=-0.1
    FODO:line=(F,O,D,O)''')
    assert_equal(1,len(mp.lines))
    fodo = mp.lines['fodo']
    assert_equal(4,len(fodo))
    assert_equal('f',fodo[0])
    assert_equal('o',fodo[1])
    assert_equal('d',fodo[2])
    assert_equal('o',fodo[3])

def test_continuation():
    mp = Mad8_parser()
    mp.parse('''q1: quadrupole,l=&
    3.14,k1=0.2''')
    assert_equal(1,len(mp.labels))
    key = 'q1'
    command = mp.labels[key]
    assert_equal('quadrupole',command.name)
    attributes = command.attributes
    assert_equal(2,len(attributes))
    assert_almost_equal(3.14,attributes['l'])
    assert_almost_equal(0.2,attributes['k1'])
