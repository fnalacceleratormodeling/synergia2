#!/usr/bin/env python

import sys
sys.path.append('..')

from nose.tools import *
from mad8_reader import Mad8_reader
from pylattice import Lattice_element, Element_adaptor_map, Lattice

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
    assert_equal(drift.get_name(),"o")
