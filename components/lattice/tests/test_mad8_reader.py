#!/usr/bin/env python

import sys
sys.path.append('..')
sys.path.append('../../foundation')

from nose.tools import *
from mad8_reader import Mad8_reader
from pylattice import Lattice_element, Element_adaptor_map, Lattice
from pyfoundation import Four_momentum, Reference_particle, constants

def test_construct():
    reader = Mad8_reader()

def test_construct2():
    adaptor_map = Element_adaptor_map()
    reader = Mad8_reader(adaptor_map)

def test_get_element_extractor_map():
    reader = Mad8_reader()
    adaptor_map = reader.get_element_adaptor_map()

def test_parse_string():
    reader = Mad8_reader()
    reader.parse_string("d:drift,l=1.2")

def test_parse():
    reader = Mad8_reader()
    reader.parse("fodo.lat")

def test_get_lines1():
    reader = Mad8_reader()
    caught = False
    try:
        lines = reader.get_lines()
    except RuntimeError, e:
        caught = True
    assert caught

def test_get_lines2():
    reader = Mad8_reader()
    reader.parse("fodo.lat")
    lines = reader.get_lines()
    lines.sort()
    expecteds = ["justf", "justo", "justd", "fodo", "model"]
    expecteds.sort()
    for line, expected in zip(lines, expecteds):
        assert_equal(line, expected)

def test_get_lattice_element1():
    reader = Mad8_reader()
    caught = False
    try:
        drift = reader.get_lattice_element("o")
    except RuntimeError, e:
        caught = True
    assert caught

def test_get_lattice_element2():
    reader = Mad8_reader()
    reader.parse("fodo.lat")
    drift = reader.get_lattice_element("o")
    assert_equal(drift.get_name(), "o")

def test_get_lattice1():
    reader = Mad8_reader()
    caught = False
    try:
        reader.get_lattice("fodo")
    except RuntimeError, e:
        caught = True
    assert caught

def assert_element_names(elements, names):
    assert_equal(len(elements), len(names))
    for element, name in zip(elements, names):
        assert_equal(element.get_name(), name)

def test_get_lattice2():
    reader = Mad8_reader()
    lattice = reader.get_lattice("fodo", "fodo.lat")
    elements = lattice.get_elements()
    assert_element_names(elements, ['f', 'o', 'd', 'o'])

def test_get_lattice_with_multiplier():
    reader = Mad8_reader()
    reader.parse_string('''
    a: drift, l=1.0
    b: drift, l=1.0
    c: drift, l=1.0
    one: line = (a,b)
    two: line = (c)
    three: line = (3*one,two)
    ''')
    lattice = reader.get_lattice("three")
    elements = lattice.get_elements()
    assert_element_names(elements, ['a', 'b', 'a', 'b', 'a', 'b', 'c'])

def test_get_lattice_with_multiplier2():
    reader = Mad8_reader()
    reader.parse_string('''
    a: drift, l=1.0
    b: drift, l=1.0
    c: drift, l=1.0
    two: line = (c)
    three: line = (3*(a,b),two)
    ''')
    lattice = reader.get_lattice("three")
    elements = lattice.get_elements()
    assert_element_names(elements, ['a', 'b', 'a', 'b', 'a', 'b', 'c'])

def test_get_lattice_with_subgroup():
    reader = Mad8_reader()
    reader.parse_string('''
    a: drift, l=1.0
    b: drift, l=1.0
    c: drift, l=1.0
    test: line = ((a,b),c)
    ''')
    lattice = reader.get_lattice("test")
    elements = lattice.get_elements()
    assert_element_names(elements, ['a', 'b', 'c'])

def test_get_lattice_with_neg():
    reader = Mad8_reader()
    reader.parse_string('''
    a: drift, l=1.0
    b: drift, l=1.0
    forw: line = (a,b)
    back: line = (-forw)
    ''')
    lattice = reader.get_lattice("back")
    elements = lattice.get_elements()
    assert_element_names(elements, ['b', 'a'])

# The example from Section 4.1.3 of the Mad8 User Guide
def test_get_lattice_ug_413():
    reader = Mad8_reader()
    reader.parse_string('''
    a: drift, l=1.0
    b: drift, l=1.0
    c: drift, l=1.0
    d: drift, l=1.0
    e: drift, l=1.0
    f: drift, l=1.0
    g: drift, l=1.0
    h: drift, l=1.0
    R: LINE=(G,H)
    S: LINE=(C,R,D)
    T: LINE=(2*S,2*(E,F),-S,-(A,B))
    ''')
    lattice = reader.get_lattice("t")
    elements = lattice.get_elements()
    assert_element_names(elements,
                        ['c', 'g', 'h', 'd', 'c', 'g', 'h', 'd', 'e', 'f', 'e', 'f', 'd', 'h', 'g', 'c', 'b', 'a'])

def test_quad_tilt():
    reader = Mad8_reader()
    reader.parse_string('q: quad, l=1.0, tilt, k1=3.14')
    element = reader.get_lattice_element('q')
    assert_equal(element.get_string_attribute('tilt'), '')
    assert(not element.has_double_attribute('tilt'))

def test_reference_particle_none():
    reader = Mad8_reader()
    reader.parse_string('''
    a: drift, l=1.0
    b: drift, l=1.0
    one: line = (a,b)
    ''')
    lattice = reader.get_lattice("one")
    assert(not lattice.has_reference_particle())

precision = 1.0e-12
def test_reference_particle():
    reader = Mad8_reader()
    lattice = reader.get_lattice("fodo", "fodo.lat")
    assert(lattice.has_reference_particle())
    reference_particle = lattice.get_reference_particle()
    four_momentum = Four_momentum(constants.mp)
    four_momentum.set_total_energy(1.5)
    expected_rp = Reference_particle(four_momentum)
    assert(reference_particle.equal(expected_rp, precision))
