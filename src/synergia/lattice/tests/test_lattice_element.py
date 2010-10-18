#!/usr/bin/env python

import sys
sys.path.append('..')
sys.path.append('../../convertors')

from pylattice import Lattice_element
from nose.tools import *

name = "foo"
type = "bar"
attr = "baz"
dblval = 2.71828
strval = "qux"

def test_construct():
    lattice_element = Lattice_element(type, name)

def test_get_type():
    lattice_element = Lattice_element(type, name)
    assert_equal(type, lattice_element.get_type())

def test_get_name():
    lattice_element = Lattice_element(type, name)
    assert_equal(name, lattice_element.get_name())

def test_has_double_attribute():
    lattice_element = Lattice_element(type, name)
    assert_equal(lattice_element.has_double_attribute(attr), False)

def test_get_set_double_attribute():
    lattice_element = Lattice_element(type, name)
    lattice_element.set_double_attribute(attr, dblval)
    assert_equal(lattice_element.has_double_attribute(attr), True)
    assert_almost_equal(lattice_element.get_double_attribute(attr), dblval)

def test_has_string_attribute():
    lattice_element = Lattice_element(type, name)
    assert_equal(lattice_element.has_string_attribute(attr), False)

def test_get_set_string_attribute():
    lattice_element = Lattice_element(type, name)
    lattice_element.set_string_attribute(attr, strval)
    assert_equal(lattice_element.has_string_attribute(attr), True)
    assert_equal(lattice_element.get_string_attribute(attr), strval)

def test_add_ancestor():
    lattice_element = Lattice_element(type, name)

    lattice_element.add_ancestor("pa")
    lattice_element.add_ancestor("grandpa")
    ancestors = lattice_element.get_ancestors()
    assert_equal(ancestors[0], "pa")
    assert_equal(ancestors[1], "grandpa")


