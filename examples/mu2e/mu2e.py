#!/usr/bin/env python
import sys
import time
import synergia
import numpy as np
from synergia.optics.one_turn_map import linear_one_turn_map
from mu2e_options import opts
import mpi4py.MPI as MPI
from basic_toolkit import *
from mxyzptlk import *
from beamline import *
from physics_toolkit import *

#####################################

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0] + csmap[1,1])
    asinmu = 0.5 * (csmap[0,0] - csmap[1,1])

    if abs(cosmu) > 1.0:
        print "error, map is unstable"
    mu = np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

#######################################################

real_particles = opts.real_particles
gridx = opts.gridx
gridy = opts.gridy
gridz = opts.gridz
partpercell = opts.partpercell

macro_particles = gridx * gridy * gridz * partpercell

seed = 4
grid = [gridx, gridy, gridz]
num_steps = opts.num_steps
num_turns = opts.num_turns
map_order = opts.map_order

radius = opts.radius
emitx = opts.norm_emit
emity = opts.norm_emit
stdz = opts.stdz
bunchlen_sec = opts.bunchlen * 1e-9
rf_voltage = opts.rf_voltage

x_offset = opts.x_offset
y_offset = opts.y_offset
z_offset = opts.z_offset

verbose = opts.verbose

myrank = MPI.COMM_WORLD.rank
if myrank == 0:
    print "==== Run Summary ===="
    print "Macro particles: ", macro_particles
    print "Real particles: ", real_particles
    print "num_turns: ", num_turns
    print "num steps/turn: ", num_steps
    print "grid: ", grid
    print "map order: ", map_order
    print "(normalized) transverse emittance: ", opts.norm_emit
    print "stdz: ", stdz
    print "bunchlen (in sec): ", bunchlen_sec
    print "Radius cut: ", radius
    print "RF Voltage: ", rf_voltage
    print "X offset: ", x_offset
    print "Y offset: ", y_offset
    print "Z_offset: ", z_offset
    print "Space Charge is ", 
    if opts.spacecharge:
        print "ON"
    else:
        print "not on"
    print "Generating diagnostics: ",
    if opts.turn_full2:
        print "turn_full2 ",
    if opts.turn_particles:
        print "turn_particles ",
    if opts.turn_tracks:
        print "turn_tracks:", opts.turn_tracks,
    print

harmno = 4

lattice = synergia.lattice.Mad8_reader().get_lattice("debunch", 
                "Debunch_modified_rf.lat")
orig_lattice_length = lattice.get_length()
if myrank == 0:
    print "original lattice length: ", orig_lattice_length

orig_elements = lattice.get_elements()
#lattice = synergia.lattice.Lattice()
#lattice.set_reference_particle(lattice.get_reference_particle())

# with the aperture, all the particles are immediately eliminated
if radius > 0.0:
    for elem in lattice.get_elements():
        elem.set_double_attribute("aperture_radius", radius)

lattice_length = lattice.get_length()

reference_particle = lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
brho = reference_particle.get_momentum() / (synergia.foundation.pconstants.c / 1e9)
bunchlen_m = bunchlen_sec * beta * synergia.foundation.pconstants.c

emitx /= (beta * gamma) * 6.0
emity /= (beta * gamma) * 6.0

if myrank == 0:
    print "lattice length: ", lattice_length
    print "brho: ", brho
    print "energy: ", energy
    print "beta: ", beta
    print "gamma: ", gamma
    print "bunchlen (in m): ", bunchlen_m

###############################################################################
#   print beamline element name, type, strength, and length
###############################################################################
#if myrank == 0:
#    elm_number = 0
#    arc_length = 0
#    for element in orig_elements:
#        elm_number +=1
#        elm_name = element.get_name()
#        elm_type = element.get_type()
#        elm_length = element.get_length()
#        #if (element.has_double_attribute(elm_name)):
#        if (elm_type == "quadrupole"):
#            elm_strength = element.get_double_attribute("k1") * brho
#        elif (elm_type == "sextupole"):
#            elm_strength = element.get_double_attribute("k2") * brho / 2.0
#        elif (elm_type == "thinSextupole"):
#            elm_strength = element.get_double_attribute("k2l") * brho / 2.0
#        elif (elm_type == "multipole"):
#            elm_strength = element.get_double_attribute("k2l") * brho / 2.0
#        elif (elm_type == "sbend"):
#            elm_strength = element.get_double_attribute("angle") * brho / elm_length
#        else:
#            elm_strength = 0
#        arc_length += elm_length
#        print elm_number, elm_name, elm_type, elm_strength, elm_length
#        if elm_name == "e_septum":
#            print "found e_septum", element.get_type()
#        if elm_name == "lambertson":
#            print "found lambertson", element.get_type()
###############################################################################



##############################################################################
#   rf cavity is not implemented yet (FIXME)
##############################################################################
# set rf cavity frequency
# harmno * beta * c/ring_length
freq = harmno * beta * synergia.foundation.pconstants.c/lattice_length
if myrank == 0:
    print "RF frequency: ", freq
if myrank == 0:
    print "Begin setting RF voltage..."
for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", rf_voltage)
        elem.set_double_attribute("freq", freq)
if myrank == 0:
    print "Finish setting RF voltage..."
###############################################################################

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, map_order)

chef_lattice = lattice_simulator.get_chef_lattice()
#print "dir(chef_lattice): ", dir(chef_lattice)

bml = chef_lattice.get_beamline()
#if myrank == 0:
#    print synergia.lattice.print_chef_beamline(bml)

proton = Proton(energy)
proton.setStateToZero()
#if myrank == 0:
#    print "original proton state: ", proton.State()
#bml.propagate(proton)
#if myrank == 0:
#    print "propagated proton: ", proton.State()

###############################################################################
#  CHEF One Turn Map
###############################################################################
JetParticle.createStandardEnvironments(map_order)
jpr = JetProton(energy)
#print "after create jpr"
bml.propagate(jpr)
chefmap = jpr.State()
#print "chefmap: ", chefmap
chefoneturnmap = chefmap.jacobian()
#print "chefoneturnmap: ", chefoneturnmap
#sys.exit(10)
###############################################################################


#chef_part = synergia.lattice.reference_particle_to_chef_particle(reference_particle)
#print "dir(chef_part): ", dir(chef_part)
#sys.exit(10)
#print "orig chef_part state: ", chef_part.State()
#chef_lattice.propagate(chef_part)
#print "propagated chef_part state: ", chef_part.State()

#sys.exit(10)
###############################################################################

map = linear_one_turn_map(lattice_simulator)
if myrank == 0:
    print "one turn map from synergia2.5 infrastructure"
    print np.array2string(map, max_line_width=200, precision=3)
    print "det(M) = ", np.linalg.det(map)


#[l, v] = np.linalg.eig(map)
#print "l: ", l
#print "v: ", v

#if myrank == 0:
#    print "eigenvalues: "
#    for z in l:
#        print "|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi)

[ax, bx, qx] = map2twiss(map[0:2, 0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6, 4:6])
dpop = bunchlen_m/bz

#alpha, beta = synergia.optics.get_alpha_beta(map)

if myrank == 0:
    print "Lattice parameters (assuming uncoupled map)"
    print "alpha_x: ", ax, " alpha_y: ", ay
    print "beta_x: ", bx, " beta_y: ", by
    print "q_x: ", qx, " q_y: ", qy
    print "beta_z: ", bz
    print "delta p/p: ", dpop

###############################################################################
#   Testing lattice with reference particles
###############################################################################
#pr1 = Proton(energy)
#pr1p = Proton(energy)
#pr2 = Proton(energy)
#pr2p = Proton(energy)
#pr3 = Proton(energy)
#pr3p = Proton(energy)
#pr1.set_x(1.0e-4)
#pr1p.set_npx(1.0e-4/bx)
#pr2.set_y(1.0e-4)
#pr2p.set_npy(1.0e-4/by)
#
#if myrank == 0:
#    print "Before propagating"
#    print "propagated pr1: ", pr1.State()
#    print "propagated pr1p: ", pr1p.State()
#    print "propagated pr2: ", pr2.State()
#    print "propagated pr2p: ", pr2p.State()
#bml.propagate(pr1)
#bml.propagate(pr1p)
#bml.propagate(pr2)
#bml.propagate(pr2p)
#if myrank == 0:
#    print "After propagating"
#    print "propagated pr1: ", pr1.State()
#    print "propagated pr1p: ", pr1p.State()
#    print "propagated pr2: ", pr2.State()
#    print "propagated pr2p: ", pr2p.State()
###############################################################################

emittance = [emitx, emity]
alpha = [ax, ay]
beta = [bx, by]
width, r = synergia.optics.match_transverse_twiss_emittance(emittance, alpha, beta)
if myrank == 0:
    print "Match transverse twiss emittance"
    print "emitx = ", emitx, "emity = ", emity
    print "xwidth = ", np.sqrt(emitx * bx), "ywidth = ", np.sqrt(emity * by)

if myrank == 0:
    print "Begin generating bunch..."

#bunch = synergia.optics.generate_matched_bunch_transverse(lattice_simulator, 
#                emitx, emity, stdz, dpop, real_particles, macro_particles,
#                seed=seed)

bunch = synergia.optics.generate_matched_bunch(lattice_simulator, 
                np.sqrt(emitx * bx), np.sqrt(emity * by), bunchlen_m, 
                real_particles, macro_particles, seed=seed)

# apply offset to bunch
particles = bunch.get_local_particles()

particles[:,0] = particles[:,0] + x_offset
particles[:,2] = particles[:,2] + y_offset
particles[:,4] = particles[:,4] + z_offset

if myrank == 0:
    print "Finish generating bunch"

if myrank == 0:
    print "expected stdx: ", np.sqrt(emitx * bx), " generated (on rank 0): ", np.std(particles[:,0])
    print "expected stdy: ", np.sqrt(emity * by), " generated (on rank 0): ", np.std(particles[:,2])
    print "expected stdz: ", bunchlen_m, " generated (on rank 0): ", np.std(particles[:,4])


###############################################################################
#   Collective operator
###############################################################################
solver = opts.spacecharge
if solver == "3d" or solver == "3D":
    if myrank == 0:
        print "using 3D Open Hockney space charge solver"
    space_charge = synergia.collective.Space_charge_3d_open_hockney(
                    bunch.get_comm(), grid)
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, 
                    space_charge, num_steps)
elif solver == "2d" or solver == "2D":
    if myrank == 0:
        print "using 2D Open Hockney space charge solver"
    space_charge = synergia.collective.Space_charge_2d_open_hockney(
                    bunch.get_comm(), grid)
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, 
                    space_charge, num_steps)
else:
    no_op = synergia.simulation.Dummy_collective_operator("stub")
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, 
                    no_op, num_steps)
###############################################################################

###############################################################################
#   Diagnostics
###############################################################################
multi_diagnostics_step = synergia.bunch.Multi_diagnostics()
# not doing step diagnostics
#multi_diagnostics_step.append(synergia.bunch.Diagnostics_full2(bunch, "mu2e_step_full2.h5"))

multi_diagnostics_turn = synergia.bunch.Multi_diagnostics()

if opts.turn_full2:
    multi_diagnostics_turn.append(synergia.bunch.Diagnostics_full2(bunch, 
                    "mu2e_full2.h5"))

if opts.turn_particles:
    multi_diagnostics_turn.append(synergia.bunch.Diagnostics_particles(bunch, 
                    "mu2e_particles.h5"))
###############################################################################

###############################################################################
#   Propagate and track particles
###############################################################################
for part in range(0, opts.turn_tracks):
    multi_diagnostics_turn.append(synergia.bunch.Diagnostics_track(bunch, 
                    "turn_track_%02d.h5"%part,part))

if myrank == 0:
    print "Begin propagate..."

t0 = time.time()
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, num_turns, multi_diagnostics_step,
        multi_diagnostics_turn, verbose)
t1 = time.time()

if myrank == 0:
    print "propagate time =", t1 - t0

if myrank == 0:
    print "Finish propagate"
