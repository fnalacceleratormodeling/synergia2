#!/usr/bin/env python
import sys
import time
import synergia
import numpy
import re
import math
from synergia.optics.one_turn_map import linear_one_turn_map
from protoplasma_options import opts
import mpi4py.MPI as MPI
from basic_toolkit import *
from mxyzptlk import *
from beamline import *
from physics_toolkit import *

###############################################################################
###############################################################################
#                                                                             #
#   ProtoPlsma Simulation in Tevatron                                         #
#                                                                             #
###############################################################################
###############################################################################
t0_start = time.time()
gridx = opts.gridx
gridy = opts.gridy
gridz = opts.gridz
partpercell = opts.partpercell
real_particles = opts.real_particles
if opts.macro_particles:
    macro_particles = opts.macro_particles
else:
    macro_particles = gridx * gridy * gridz * partpercell

seed = 4
grid = [gridx, gridy, gridz]
num_steps = opts.num_steps
num_turns = opts.num_turns
map_order = opts.map_order

radius = opts.radius
sigx = opts.sigx
sigy = opts.sigy
bunchlen_sec = opts.bunchlen * 1e-9

x_offset = opts.x_offset
y_offset = opts.y_offset
z_offset = opts.z_offset

solver = opts.spacecharge
verbose = opts.verbose

myrank = MPI.COMM_WORLD.rank
if myrank == 0:
    print
    print "Run Summary"
    print "    Macro particles                  :", macro_particles
    print "    Real particles                   :", real_particles
    print "    num_turns                        :", num_turns
    print "    num steps/turn                   :", num_steps
    print "    grid                             :", grid
    print "    map order                        :", map_order
    print "    bunch length                     :", bunchlen_sec * 1e9, "nsec"
    print "    Radius aperture                  :", radius, "m"
    print "    X offset                         :", x_offset
    print "    Y offset                         :", y_offset
    print "    Z_offset                         :", z_offset
    print "    Space Charge                     :",
    if solver == "2d" or solver == "2D" or solver == "3d" or solver == "3D":
        print "ON"
    else:
        print "OFF"
    print "    Generating diagnostics           :",
    if opts.turn_full2:
        print "turn_full2",
    if opts.turn_particles:
        print "turn_particles",
    if opts.turn_tracks:
        print "turn_tracks:", opts.turn_tracks,
    print

synergia_lattice = synergia.lattice.Mad8_reader().get_lattice(
                "protoplasma", "Tevatron-E.lat")
synergia_elements = synergia_lattice.get_elements()

tevatron_lattice = synergia.lattice.Mad8_reader().get_lattice(
                "tevatron", "Tevatron-E.lat")

lattice_length = synergia_lattice.get_length()

reference_particle = synergia_lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
momentum = reference_particle.get_momentum()
brho = momentum / (synergia.foundation.pconstants.c / 1e9)
bunchlen_m = bunchlen_sec * beta * synergia.foundation.pconstants.c

if myrank == 0:
    print "    lattice length                   :", lattice_length, "m"
    print "    brho                             :", brho, "T-m"
    print "    momentum                         :", momentum, "GeV/c"
    print "    energy                           :", energy, "GeV"
    print "    beta                             :", beta
    print "    gamma                            :", gamma
    print "    bunch length (in second)         :", bunchlen_sec * 1e9, "nsec"
    print "    bunch length (in meter)          :", bunchlen_m, "m"

###############################################################################
#   Set lattice_simulator
###############################################################################
if myrank == 0:
    print
    print "Set lattice_simulator"

lattice_simulator = synergia.simulation.Lattice_simulator(synergia_lattice, 
                map_order)
tevatron_lattice_simulator = synergia.simulation.Lattice_simulator(
                tevatron_lattice, map_order)

###############################################################################
#   Set/Store initial element strengths
###############################################################################
'''
if myrank == 0:
    print
    print "Set initial element strengths"
    print
    print "....Set initial settings for the Synergia lattice"
index = 0
for element in synergia_elements:
    if element.get_type() == "quadrupole":
        element.set_double_attribute("k1", initial_k1[index])
        index += 1
    if element.get_type() == "multipole":
        element.set_double_attribute("k2l", 0.0)
'''
lattice_simulator.update()

###############################################################################
#   Lattice functions
###############################################################################
'''
index = 0
for element in synergia_elements:
    lattice_functions = lattice_simulator.get_lattice_functions(element)
    beta_x = lattice_functions.beta_x
    alpha_x = lattice_functions.alpha_x
    D_x = lattice_functions.D_x
    if myrank == 0: 
        print
        print "        index                        :", index
        print "        type                         :", 
        print element.get_type()
        print "        name                         :", 
        print element.get_name()
        print "        beta_x                       :", beta_x
        print "        alpha_x                      :", alpha_x
        print "        D_x                          :", D_x
    index += 1
'''

#   Output files for plotting
twiss_file = ("twiss.txt")
twiss_log = open(twiss_file, "w")
index = 0
length = 0.0
for element in synergia_elements:
    lattice_functions = lattice_simulator.get_lattice_functions(element)
    type = element.get_type()
    name = element.get_name()
    length += element.get_length()
    beta_x = lattice_functions.beta_x
    alpha_x = lattice_functions.alpha_x
    D_x = lattice_functions.D_x
    if myrank == 0:
        twiss_log.write("%3d %12s %10s %10.5f %10.5f %10.5f %10.5f\n" % (
                        index, type, name, length, beta_x, alpha_x, D_x))
        twiss_log.flush()
    index += 1
twiss_log.close()
###############################################################################
#   Generating Bunch
###############################################################################
if myrank == 0:
    print
    print "Begin generating bunch..."
bunch = synergia.optics.generate_matched_bunch(lattice_simulator,
        sigx, sigy, bunchlen_m, real_particles, macro_particles, seed=seed)
#bunch = synergia.optics.generate_matched_bunch(tevatron_lattice_simulator,
#            sigx, sigy, bunchlen_m, real_particles, macro_particles, seed=seed)

# apply offset to bunch
particles = bunch.get_local_particles()

particles[:,0] = particles[:,0] + x_offset
particles[:,2] = particles[:,2] + y_offset
particles[:,4] = particles[:,4] + z_offset

###############################################################################
#   Collective operator
###############################################################################
if myrank == 0:
    print
    print "Set collective operator"
if solver == "3d" or solver == "3D":
    if myrank == 0:
        print "    using 3D Open Hockney space charge solver"
    space_charge = synergia.collective.Space_charge_3d_open_hockney(
                    bunch.get_comm(), grid, False)
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator,
                    space_charge, num_steps)
elif solver == "2d" or solver == "2D":
    if myrank == 0:
        print "    using 2D Open Hockney space charge solver"
    space_charge = synergia.collective.Space_charge_2d_open_hockney(
                    bunch.get_comm(), grid)
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator,
                    space_charge, num_steps)
else:
    if myrank == 0:
        print "    no collective operator is selected"
    no_op = synergia.simulation.Dummy_collective_operator("stub")
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator,
                    no_op, num_steps)

###############################################################################
#   Diagnostics
###############################################################################
diagnostics_actions = synergia.simulation.Diagnostics_actions()
#diagnostics per step
for part in range(0, opts.step_tracks):
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_track(bunch,
                            "tevatron_step_track_%02d.h5" % part, part))
if opts.step_full2:
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_full2(bunch,
                            "tevatron_step_full2.h5"))
if opts.step_particles:
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_particles(bunch,
                            "tevatron_step_particles.h5"))
# diagnostics per turn
for part in range(0, opts.turn_tracks):
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_track(bunch,
                            "tevatron_turn_track_%02d.h5" % part, part))
if opts.turn_full2:
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_full2(bunch,
                            "tevatron_turn_full2.h5"))
if opts.turn_particles:
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_particles(bunch,
                            "tevatron_turn_particles.h5"))

###############################################################################
#   Propagate and track particles
###############################################################################
if myrank == 0:
    print
    print "Begin propagate..."

t0 = time.time()
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, num_turns, diagnostics_actions,  verbose)
t1 = time.time()
