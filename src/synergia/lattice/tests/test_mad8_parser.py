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
    except ParseException as e:
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

def test_variable_assignment():
    mp = Mad8_parser()
    # 7/32 has an exact representation in floating point
    mp.parse('x=0.21875E+02')
    assert_equal(21.875,mp.variables['x'])
    mp.parse('y=.21875E+02')
    assert_equal(21.875,mp.variables['y'])

def test_caps_variable_assignment():
    mp = Mad8_parser()
    mp.parse('X=1;Y=X')
    assert_equal(1,mp.variables['y'])

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

def test_subscripted_ident():
    mp = Mad8_parser()
    mp.parse('foo: bar, a=1;x=foo[a]')
    assert_equal(1,mp.variables['x'])

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

def test_continuation2():
    mp = Mad8_parser()
    mp.parse('''q2: quadrupole,l=& )junk)
    3.14,k1=0.2''')
    assert_equal(1,len(mp.labels))
    key = 'q2'
    command = mp.labels[key]
    assert_equal('quadrupole',command.name)
    attributes = command.attributes
    assert_equal(2,len(attributes))
    assert_almost_equal(3.14,attributes['l'])
    assert_almost_equal(0.2,attributes['k1'])

def test_continuation3():
    mp = Mad8_parser()
    mp.parse('''q3: quadrupole,l=3.14 &
    k1=0.2''')
    assert_equal(1,len(mp.labels))
    key = 'q3'
    command = mp.labels[key]
    assert_equal('quadrupole',command.name)
    attributes = command.attributes
    assert_equal(2,len(attributes))
    assert_almost_equal(3.14,attributes['l'])
    assert_almost_equal(0.2,attributes['k1'])

def test_continuation4():
    mp = Mad8_parser()
    mp.parse('''twiss, save,   betx=28.871,alfx=-0.069,mux=0.0,dx=2.682,dpx=-0.073 &
               bety= 5.264,alfy=-0.006,muy=0.0,dy=0.0,dpy=0.0''')

def test_line_multiplier():
    mp = Mad8_parser()
    mp.parse('fodo8: line=(8*(f,o,d,o))')
    print(mp.lines)
    assert_equal(1,len(mp.lines))
    assert_equal(2,len(mp.lines['fodo8']))
    assert_equal('8*',mp.lines['fodo8'][0])
    assert_equal(4,len(mp.lines['fodo8'][1]))
    assert_equal('f',mp.lines['fodo8'][1][0])
    assert_equal('o',mp.lines['fodo8'][1][1])
    assert_equal('d',mp.lines['fodo8'][1][2])
    assert_equal('o',mp.lines['fodo8'][1][3])

def test_line_unary_minus():
    mp = Mad8_parser()
    mp.parse('odof: line=(-(f,o,d,o))')
    print(mp.lines)
    assert_equal(1,len(mp.lines))
    assert_equal(2,len(mp.lines['odof']))
    assert_equal('-',mp.lines['odof'][0])
    assert_equal(4,len(mp.lines['odof'][1]))
    assert_equal('f',mp.lines['odof'][1][0])
    assert_equal('o',mp.lines['odof'][1][1])
    assert_equal('d',mp.lines['odof'][1][2])
    assert_equal('o',mp.lines['odof'][1][3])

def test_range1():
    mp = Mad8_parser()
    mp.parse('print, #s/#e')

def test_range2():
    mp = Mad8_parser()
    mp.parse('select,optics #s/#e')

def test_multi_value():
    mp = Mad8_parser()
    mp.parse('''plot, vmin=0,0''')
