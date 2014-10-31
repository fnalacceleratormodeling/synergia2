# !/usr/bin/env synergia
import synergia
import numpy as np
import sys

from synergia.foundation import Four_momentum, Reference_particle, pconstants
from synergia.lattice import Mad8_reader, Lattice
from synergia.bunch import Bunch, Diagnostics_basic, Diagnostics_track, Diagnostics_particles
from synergia.simulation import Independent_stepper_elements, Bunch_simulator, \
    Propagator

# Get options from separate options file
from fodo_set_beam_by_hand_options import opts

# Define a lattice
#     Read the lattice named "fodo" from the Mad8 file "fodo.lat"
lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")

refpart = lattice.get_reference_particle()
energy = refpart.get_total_energy()
momentum = refpart.get_momentum()
gamma = refpart.get_gamma()
beta = refpart.get_beta()
print "Beamline energy: ", energy
print "Beamline momentum: ", momentum
print "Beamline particle gamma: ", gamma
print "Beamline particle beta: ", beta

# Define a set of simulation steps
stepper = Independent_stepper_elements(lattice, opts.map_order, 
                                       opts.steps_per_element)

# Make a and populate bunch by hand

# first create the bunch

# I need a communicator so I can make a proper  multi-processor distribution
commxx = synergia.utils.Commxx()

bunch = synergia.bunch.Bunch(
    lattice.get_reference_particle(),
    opts.macro_particles, opts.real_particles, commxx)

# Each processor gets its own array of particles
local_particles = bunch.get_local_particles()

# my current processor has local_num particles
local_num = bunch.get_local_num()

# Now I can fill my local particles however I want by filling
# local_particles[partnum, coord] coord=0:5 {x, xp, y, yp, cdt, dpop]
# and partnum is 0:local_num-1

# get lattice functions to help set the distribution
last_lattice_element = lattice.get_elements()[-1]

lf = stepper.get_lattice_simulator().get_lattice_functions(last_lattice_element)

beta_x = lf.beta_x
alpha_x = lf.alpha_x
beta_y = lf.beta_y
alpha_y = lf.alpha_y
print "beta x: ", beta_x
print "alpha x: ", alpha_x
print "beta y: ", beta_y
print "alpha y: ", alpha_y

# generate a hollow transverse distribution
# x and (alpha*x+beta*xp) are the normal components of the distribution

dist = synergia.foundation.Random_distribution(opts.seed, commxx)

rx = opts.x_std*np.sqrt(2.0)

# fill the x and xp components being the two components of the
# normal form x, and alpha*x+beta*xp
theta1 = np.zeros(local_num, 'd')
dist.fill_uniform(theta1, 0.0, 2.0*np.pi)
local_particles[:,0] = rx * np.cos(theta1)
local_particles[:,1] = (rx*np.sin(theta1) - alpha_x*local_particles[:,0])/beta_x

print "x std: ", local_particles[:,0].std()
print "xp std: ", local_particles[:,1].std()

# do the same for y
ry = opts.y_std*np.sqrt(2.0)

theta2 = np.zeros(local_num, 'd')
dist.fill_uniform(theta2, 0.0, 2.0*np.pi)

local_particles[:,2] = ry * np.cos(theta2)
local_particles[:,3] = (ry * np.sin(theta2) - alpha_y*local_particles[:,2])/beta_y
print "y_std: ", local_particles[:,2].std()
print "yp std: ", local_particles[:,3].std()

# longitudinal is a straight gaussian

zdist = np.zeros(local_num, 'd')
dist.fill_unit_gaussian(zdist)
local_particles[:, 4] = zdist * opts.z_std/beta
print "c*dt std: ", local_particles[:,4].std()

# delta p/p is a straight gaussian
dpopdist = np.zeros(local_num, 'd')
dist.fill_unit_gaussian(dpopdist)
local_particles[:, 5] = dpopdist * opts.dpop
print "dpop std: ", local_particles[:,5].std()

# Define a bunch simulator
bunch_simulator = Bunch_simulator(bunch)

# Define a set of bunch diagnostics
#     Apply basic diagnostics every step
diagnostics = Diagnostics_basic("diagnostics.h5")
bunch_simulator.add_per_step(diagnostics)

# Optionally, save track of a single particle
if opts.track_particle:
    track_diagnostics = Diagnostics_track("track.h5", opts.track_particle)
    bunch_simulator.add_per_step(track_diagnostics)

if opts.save_particles:
    particles_diagnostics = Diagnostics_particles("particles.h5")
    bunch_simulator.add_per_turn(particles_diagnostics)

# Perform the simulation
propagator = Propagator(stepper)
propagator.propagate(bunch_simulator, opts.turns, opts.max_turns, 
                     opts.verbosity)
