#!/usr/bin/env python3
from __future__ import print_function

import sys
import os
import numpy
import synergia
from ramp_module import Ramp_actions
import re

from channel_options import opts

f = open('channel_template.seq', 'r')
channel_template = f.read()
f.close()

channel_template = re.sub('%skew', repr(opts.skew), channel_template)
channel_template = re.sub('%sext', repr(opts.sext), channel_template)
channel_template = re.sub('%sksext', repr(opts.sksext), channel_template)
channel_template = re.sub('%octo', repr(opts.octo), channel_template)
channel_template = re.sub('%skocto', repr(opts.skocto), channel_template)
channel_template = re.sub('%RFVolt', repr(opts.RFVolt), channel_template)

f = open('channel.seq', 'w')
f.write(channel_template)
f.close()

reader = synergia.lattice.MadX_reader()
reader.parse(channel_template)
lattice = reader.get_lattice("fodo")

#lattice = synergia.lattice.MadX_().get_lattice("fodo", "channel.seq")


for elem in lattice.get_elements():
    elem.set_string_attribute("extractor_type", "chef_propagate")

print("lattice: ", len(lattice.get_elements()), " elements, length: ", lattice.get_length())

reference_particle = lattice.get_reference_particle()

energy = reference_particle.get_total_energy()
momentum = reference_particle.get_momentum()
gamma = reference_particle.get_gamma()
beta = reference_particle.get_beta()
print("energy: ", energy)
print("momentum: ", momentum)
print("gamma: ", gamma)
print("beta: ", beta)

macro_particles=80
real_particles = 1.0e9

# order is the order in which the normal form calculation is done
order = opts.order
print("normal forms to order ", order)

stepper = synergia.simulation.Independent_stepper(lattice, order, 1)

lattice_simulator = stepper.get_lattice_simulator()
lattice_simulator.print_lattice_functions()

# synergia.simulation.xml_save_fnf(lattice_simulator.get_fast_normal_form(), "channel_fnf.xml")

map = lattice_simulator.get_linear_one_turn_map()
print("one turn map:")
print(numpy.array2string(map,max_line_width=200))
print()

[l, v] = numpy.linalg.eig(map)
print("eigenvalues: ")
for z in l:
    print("|z|: ", abs(z), " z: ", z, " tune: ", numpy.log(z).imag/(2.0*numpy.pi))

(nux, nuy) = lattice_simulator.get_both_tunes()
chrx = lattice_simulator.get_horizontal_chromaticity()
chry = lattice_simulator.get_vertical_chromaticity()

print("nux: ", nux)
print("nuy: ", nuy)
print("chrx: ", chrx)
print("chry: ", chry)

print("linear normal form check: ", lattice_simulator.check_linear_normal_form())

commxx = synergia.utils.Commxx()

bunch = synergia.bunch.Bunch(reference_particle, macro_particles, real_particles, commxx)

local_particles = bunch.get_local_particles()
for i in range(macro_particles//2):
    local_particles[i, 0:6] = 0.0
    local_particles[i, 0] = opts.offset*i
for i in range(macro_particles//2, macro_particles):
    local_particles[i, 0:6] = 0.0
    local_particles[i, 2] = opts.offset*(i-macro_particles/2)

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
bunch_simulator.add_per_step(synergia.bunch.Diagnostics_bulk_track("tracks.h5", macro_particles))

ramp_actions = Ramp_actions()

propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(0)
propagator.propagate(bunch_simulator, ramp_actions, opts.turns, opts.turns, 1)
