#!/usr/bin/env python

# see also: test_simplify2.py in the simulation module

import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
from synergia.lattice import eliminate_markers, convert_monitors, convert_magnets, \
    combine_drifts, simplify_all
from synergia.lattice import Lattice, xml_load_lattice

def test_booster():
    orig = Lattice()
    xml_load_lattice(orig, "lattices/fnal_booster.xml")
    simplified = simplify_all(orig)
    assert_almost_equal(orig.get_length(), simplified.get_length())
    print "booster went from", len(orig.get_elements()), "to",
    print len(simplified.get_elements()), "elements"

def test_debuncher():
    orig = Lattice()
    xml_load_lattice(orig, "lattices/fnal_debuncher.xml")
    simplified = simplify_all(orig)
    assert_almost_equal(orig.get_length(), simplified.get_length())
    print "debuncher went from", len(orig.get_elements()), "to",
    print len(simplified.get_elements()), "elements"

def test_main_injector():
    orig = Lattice()
    xml_load_lattice(orig, "lattices/fnal_main_injector.xml")
    simplified = simplify_all(orig)
    assert_almost_equal(orig.get_length(), simplified.get_length())
    print "main injector went from", len(orig.get_elements()), "to",
    print len(simplified.get_elements()), "elements"

if __name__ == '__main__':
    test_booster()
    test_debuncher()
    test_main_injector()

