#!/usr/bin/env python

# see also: test_simplify2.py in the simulation module
from __future__ import print_function
import sys
sys.path.append('../../..')
import local_paths

from nose.tools import *
from synergia.lattice import eliminate_markers, convert_monitors, convert_magnets, \
    combine_drifts, simplify_all
from synergia.lattice import Lattice
from synergia.utils import read_lsexpr_file
import synergia

def print_change(name, old, new):
    old_len = len(old.get_elements())
    new_len = len(new.get_elements())
    print(name, "we from", old_len, "to", new_len, "elements")

def test_booster():
    orig = Lattice(read_lsexpr_file("lattices/fnal_booster.lsx"))
    simplified = simplify_all(orig)
    assert_almost_equal(orig.get_length(), simplified.get_length())
    print_change("booster", orig, simplified)

def test_debuncher():
    orig = Lattice(read_lsexpr_file("lattices/fnal_debuncher.lsx"))
    simplified = simplify_all(orig)
    assert_almost_equal(orig.get_length(), simplified.get_length())
    print_change("debunche", orig, simplified)

def test_main_injector():
    orig = Lattice(read_lsexpr_file("lattices/fnal_main_injector.lsx"))
    simplified = simplify_all(orig)
    assert_almost_equal(orig.get_length(), simplified.get_length())
    print_change("main injector", orig, simplified)

if __name__ == '__main__':
    test_booster()
    test_debuncher()
    test_main_injector()

