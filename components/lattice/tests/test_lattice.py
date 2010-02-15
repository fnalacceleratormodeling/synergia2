#!/usr/bin/env python

import sys
sys.path.append('..')

from pylattice import Lattice_element
from nose.tools import *

name = "foo"
type = "bar"
attr = "baz"
val = 2.71828

def test_construct_lattice_element():
    lattice_element = Lattice_element(type, name)

def test_get_type():
    lattice_element = Lattice_element(type, name)
    assert_equal(type, lattice_element.get_type())

def test_get_name():
    lattice_element = Lattice_element(type, name)
    assert_equal(name, lattice_element.get_name())

def test_has_attribute():
    lattice_element = Lattice_element(type, name)
    assert_equal(lattice_element.has_double_attribute(attr), False)

def test_get_set_attribute():
    lattice_element = Lattice_element(type, name)
    lattice_element.set_double_attribute(attr, val)
    assert_equal(lattice_element.has_double_attribute(attr), True)
    assert_almost_equal(lattice_element.get_double_attribute(attr), val)
