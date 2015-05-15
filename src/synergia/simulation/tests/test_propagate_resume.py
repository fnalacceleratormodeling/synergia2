#!/usr/bin/env python

from nose.tools import *
import synergia
from synergia.utils import Hdf5_file
from synergia.simulation import Lattice_simulator, Bunch_simulator, Independent_stepper, Propagator, Resume
from synergia.bunch import Bunch, Diagnostics_particles, Diagnostics_full2, Diagnostics_bulk_track
from synergia.optics import generate_matched_bunch_transverse

import math

name = "fodo"
charge = 1
mass = synergia.foundation.pconstants.mp
p = 1.5
total_energy = math.sqrt(p**2 + mass**2)
tolerance = 1.0e-12
map_order = 1
n_cells = 1

focus = 7;
sepn = 10;
quad_length = 0.2;
quad_strength = 1.0 / (focus * quad_length);
pct = 0.4;
drift_length = sepn - quad_length

map_order = 1
macro_particles = 8192
real_particles = 1.0e11
emit = 1.0e-6
rmsz = 0.1
dpop = 1.0e-4

nturns = 100
maxturns = 10

class Fixture:
    def __init__(self):
        four_momentum = synergia.foundation.Four_momentum(mass, total_energy)
        reference_particle = \
            synergia.foundation.Reference_particle(charge, four_momentum)
        self.lattice = synergia.lattice.Lattice(name)
        self.lattice.set_reference_particle(reference_particle)
        f = synergia.lattice.Lattice_element("quadrupole", "f")
        f.set_double_attribute("l", quad_length)
        f.set_double_attribute("k1", quad_strength)
        f.set_string_attribute("extractor_map", "chef_propagate")
        o = synergia.lattice.Lattice_element("drift", "o")
        o.set_double_attribute("l", drift_length)
        o.set_string_attribute("extractor_map", "chef_propagate")
        d = synergia.lattice.Lattice_element("quadrupole", "d")
        d.set_double_attribute("l", quad_length)
        d.set_double_attribute("k1", -quad_strength)
        d.set_string_attribute("extractor_map", "chef_propagate")

        for cell in range(0, n_cells):
            self.lattice.append(f)
            self.lattice.append(o)
            self.lattice.append(d)
            self.lattice.append(o)

        self.lattice_simulator = Lattice_simulator(self.lattice, map_order)

        self.bunch = generate_matched_bunch_transverse(self.lattice_simulator, emit, emit, rmsz, dpop, real_particles, macro_particles)

class Lattice_fixture:
    def __init__(self):
        n_cells = 1
        four_momentum = synergia.foundation.Four_momentum(mass, total_energy)
        reference_particle = \
            synergia.foundation.Reference_particle(charge, four_momentum)
        self.lattice = synergia.lattice.Lattice(name)
        self.lattice.set_reference_particle(reference_particle)
        f = synergia.lattice.Lattice_element("quadrupole", "f")
        f.set_double_attribute("l", quad_length)
        f.set_double_attribute("k1", quad_strength)
        f.set_string_attribute("extractor_map", "chef_propagate")
        o = synergia.lattice.Lattice_element("drift", "o")
        o.set_double_attribute("l", drift_length)
        o.set_string_attribute("extractor_map", "chef_propagate")
        d = synergia.lattice.Lattice_element("quadrupole", "d")
        d.set_double_attribute("l", quad_length)
        d.set_double_attribute("k1", -quad_strength)
        d.set_string_attribute("extractor_map", "chef_propagate")

        for cell in range(0, n_cells):
            self.lattice.append(f)
            self.lattice.append(o)
            self.lattice.append(d)
            self.lattice.append(o)

def check_full2(file, turns, num_particles):
    bad = False
    h5 = Hdf5_file(file, Hdf5_file.read_only)
    if h5.get_dims("num_particles")[0] != turns:
        print "check_full2: wrong size of num_particles, ", h5.get_dims("num_particles")[0], "!=", turns
        bad = True
    if h5.get_dims("mean")[1] != turns:
        print "check_full2: wrong size of mean array: ", h5.get_dims("mean"), "!=", (6,turns)
    if h5.read_array1i('num_particles')[-1] != num_particles:
        print "check_full2: wrong number of macro_particles, ", h5.get_dims("num_particles")[-1], " != ", num_particles
        bad = True
    return bad

def check_tracks(file, turns, particles):
    bad = False
    h5 = Hdf5_file(file, Hdf5_file.read_only)
    if h5.get_dims("track_coords")[0] != particles:
        print "check_tracks: wrong number of particles: ", h5.get_dims("track_coords")[0], "!=", particles
        bad = True
    if h5.get_dims("track_coords")[2] != turns:
        print "check_tracks: wrong number of turns: ", h5.get_dims("track_coords")[2], "!=", turns
        bad = True
    return bad

def check_particles(file, turn0, turn1, particles):
    bad = False
    for i in range(turn0, turn1):
        filename = file+"_%04d.h5"%i
        h5 = Hdf5_file(filename, Hdf5_file.read_only)
        if h5.get_dims("particles")[0] != particles:
            print "check_particles: file ", filename, " wrong number of particles ", h5.get_dims("particles")[0], "!=", particles
            bad = True
    return bad

def test_propagate():
    bad = False
    f = Fixture()
    bunch_simulator = Bunch_simulator(f.bunch)
    bunch_simulator.add_per_turn(Diagnostics_particles("turn_particles.h5"))
    bunch_simulator.add_per_turn(Diagnostics_full2("turn_full2.h5"))
    bunch_simulator.add_per_turn(Diagnostics_bulk_track("turn_tracks.h5", f.bunch.get_total_num()))

    stepper = Independent_stepper(f.lattice, map_order, 1)
    propagator = Propagator(stepper)
    propagator.propagate(bunch_simulator, nturns, maxturns)
    # check that all the diagnostics were written out
    if check_full2("turn_full2.h5", maxturns+1, macro_particles):
        bad = True
    if check_tracks("turn_tracks.h5", maxturns+1, macro_particles):
        bad = True
    if check_particles("turn_particles", 0, maxturns+1, macro_particles):
        bad = True
    assert (not bad)


def test_resume():
    resume = Resume(Propagator.default_checkpoint_dir)
    resume.propagate(0, False, 0, False, 0, False)
    # check that all the diagnostics were written out
    bad = False
    if check_full2("turn_full2.h5", 2*maxturns+1, macro_particles):
        bad = True
    if check_tracks("turn_tracks.h5", 2*maxturns+1, macro_particles):
        bad = True
    if check_particles("turn_particles", maxturns+1, 2*maxturns+1, macro_particles):
        bad = True
    assert (not bad)

if __name__ == "__main__":
    logger = synergia.utils.Logger(0)
    if opts.propagate:
        if test_propagate():
            raise RuntimeError, "test_propagate() failed"
        else:
            print >>logger, "test_propagate() passed"
    if opts.resume:
        if test_resume():
            raise RuntimeError, "test_resume() failed"
        else:
            print >>logger, "test_resume() passed"
