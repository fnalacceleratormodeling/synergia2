#!/usr/bin/env python

import sys
sys.path.append('../../..')
import local_paths

from synergia.lattice import Lattice_element, Lattice
from nose.tools import *

name = "foo"
type = "bar"
attr = "baz"
attr2 = "qux"
dblval = 2.71828
dblval2 = 3.1415
strval = "quux"
strval2 = "corge"

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

def test_get_double_attributes():
    lattice_element = Lattice_element(type, name)
    lattice_element.set_double_attribute(attr, dblval)
    lattice_element.set_double_attribute(attr2, dblval2)
    double_attributes = lattice_element.get_double_attributes()
    assert_equal(len(double_attributes), 2)
    assert_almost_equal(double_attributes[attr], dblval)
    assert_almost_equal(double_attributes[attr2], dblval2)

def test_has_string_attribute():
    lattice_element = Lattice_element(type, name)
    assert_equal(lattice_element.has_string_attribute(attr), False)

def test_get_set_string_attribute():
    lattice_element = Lattice_element(type, name)
    lattice_element.set_string_attribute(attr, strval)
    assert_equal(lattice_element.has_string_attribute(attr), True)
    assert_equal(lattice_element.get_string_attribute(attr), strval)

def test_get_string_attributes():
    lattice_element = Lattice_element(type, name)
    lattice_element.set_string_attribute(attr, strval)
    lattice_element.set_string_attribute(attr2, strval2)
    string_attributes = lattice_element.get_string_attributes()
    assert_equal(len(string_attributes), 2)
    assert_equal(string_attributes[attr], strval)
    assert_equal(string_attributes[attr2], strval2)

def test_add_ancestor():
    lattice_element = Lattice_element(type, name)

    lattice_element.add_ancestor("pa")
    lattice_element.add_ancestor("grandpa")
    ancestors = lattice_element.get_ancestors()
    assert_equal(ancestors[0], "pa")
    assert_equal(ancestors[1], "grandpa")

def test_has_lattice():
    lattice_element = Lattice_element(type, name)
    has_one = lattice_element.has_lattice()
    assert_equal(has_one, False)

lattice_name = "foo_lattice"
def test_set_get_lattice():
    lattice_element = Lattice_element(type, name)
    lattice = Lattice(lattice_name)
    lattice_element.set_lattice(lattice)
    # This is a less-than-optimal test. The test in
    # test_lattice_element.cc is better.
    assert_equal(lattice_element.get_lattice().get_name(),
        lattice.get_name())

