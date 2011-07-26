#!/usr/bin/env python
import sys
import synergia
import numpy as np
from synergia.optics.one_turn_map import linear_one_turn_map
from mi_options import opts
import mpi4py.MPI as MPI
from basic_toolkit import *
from mxyzptlk import *
from beamline import *
from physics_toolkit import *

#####################################

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        print "error, map is unstable"
    mu =np.arccos(cosmu)

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
real_particles = real_particles
num_steps = opts.num_steps
num_turns = opts.num_turns
map_order = opts.map_order

radius = opts.radius
emitx = opts.norm_emit
emity = opts.norm_emit
stdz = opts.stdz
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
        print "turn_tracks:",opts.turn_tracks,
    print

harmno = 588

orig_lattice = synergia.lattice.Mad8_reader().get_lattice("ring_p_q605", "mi20-egs-thinrf.lat")
orig_lattice_length = orig_lattice.get_length()
if myrank == 0:
    print "original lattice length: ", orig_lattice_length

# go through lattice splitting each quadrupole and inserting the thin
# multipole object building the new lattice
if myrank == 0:
    print "Begin constructing lattice... "
orig_elements = orig_lattice.get_elements()
lattice = synergia.lattice.Lattice()
for elem in orig_elements:
    # quadrupole splitting for multipole insertion disabled
    if elem.get_type() == "xxxquadrupole":
        old_name = elem.get_name()
        old_length = elem.get_length()
        new_length = old_length/2.0
        
        # first half of split quadrupole
        new_elem1 = synergia.lattice.Lattice_element(elem.get_type(),
                                                     old_name+"_1")
        string_attrs = elem.get_string_attributes()
        double_attrs = elem.get_double_attributes()
        for k in string_attrs.keys():
            new_elem1.set_string_attribute(k, string_attrs[k])
        for k in double_attrs.keys():
            new_elem1.set_double_attribute(k, double_attrs[k])
        new_elem1.set_double_attribute("l", new_length)

        # second half of split quadrupole
        new_elem2 = synergia.lattice.Lattice_element(elem.get_type(),
                                                     old_name+"_2")
        for k in string_attrs.keys():
            new_elem2.set_string_attribute(k, string_attrs[k])
        for k in double_attrs.keys():
            new_elem2.set_double_attribute(k, double_attrs[k])
        new_elem2.set_double_attribute("l", new_length)

        # extract kl for enclosing quadrupole for normalizing the thin
        # multipole object
        kl = old_length * elem.get_double_attribute("k1")
        thinpole_elem = synergia.lattice.Lattice_element("thinpole",
                                                         old_name+"_poles")
        thinpole_elem.set_double_attribute("kl", kl)
        if not multipoles.has_key(old_name):
            if  myrank == 0:
                pass
                #print "Quadrupole ",old_name," has no entry in multipoles file"
            else:
                pass
        else:
            mcoeff = multipoles[old_name]
            for k in mcoeff.keys():
                thinpole_elem.set_double_attribute(k, mcoeff[k])

        lattice.append(new_elem1)
        lattice.append(thinpole_elem)
        lattice.append(new_elem2)

    else:
        lattice.append(elem)

lattice.set_reference_particle(orig_lattice.get_reference_particle())

# with the aperture, all the particles are immediately eliminated
if radius > 0.0:
    for elem in lattice.get_elements():
        elem.set_double_attribute("aperture_radius", radius)

if myrank == 0:
    print "Finished constructing lattice"

lattice_length = lattice.get_length()          

reference_particle = lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()

if myrank == 0:
    print "lattice length: ", lattice_length
    print "energy: ", energy
    print "beta: ", beta
    print "gamma: ", gamma

# set rf cavity frequency
# harmno * beta * c/ring_length
freq = harmno * beta * synergia.foundation.pconstants.c/lattice_length
if myrank == 0:
    print "RF frequency: ", freq

if myrank == 0:
    print "Begin setting RF voltage..."

# rf cavity voltage, is 1.0 MV total distributed over 18 cavities.  MAD8
# expects cavities voltages in  units of MV.
for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", rf_voltage)
        elem.set_double_attribute("freq", freq)

if myrank == 0:
    print "Finish setting RF voltage..."

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, map_order)

chef_lattice = lattice_simulator.get_chef_lattice()
#print "dir(chef_lattice): ", dir(chef_lattice)
bml = chef_lattice.get_beamline()
proton = Proton(energy)
proton.setStateToZero()
if myrank==0:
    print "original proton state: ", proton.State()
bml.propagate(proton)
if myrank==0:
    print "propagated proton: ", proton.State()

#jpr = JetProton(energy)
#print "after create jpr"
#bml.propagate(jpr)
#chefmap = jpr.State()
#print "chefmap: ", chefmap
#chefoneturnmap = Jacobian(chefmap)

#print "chefoneturnmap: ", chefoneturnmap

#sys.exit(10)

#chef_part = synergia.lattice.reference_particle_to_chef_particle(reference_particle)
#print "dir(chef_part): ", dir(chef_part)
#sys.exit(10)
#print "orig chef_part state: ", chef_part.State()
#chef_lattice.propagate(chef_part)
#print "propagated chef_part state: ", chef_part.State()

#sys.exit(10)

map = linear_one_turn_map(lattice_simulator)
if myrank==0:
    print "one turn map from synergia2.5 infrastructure"
    print np.array2string(map, max_line_width=200)


[l, v] = np.linalg.eig(map)

#print "l: ", l
#print "v: ", v

if myrank==0:
    print "eigenvalues: "
    for z in l:
        print "|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi)

[ax, bx, qx] = map2twiss(map[0:2,0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6,4:6])
dpop = stdz/bz

if myrank == 0:
    print "Lattice parameters (assuming uncoupled map)"
    print "alpha_x: ", ax, " alpha_y: ", ay
    print "beta_x: ", bx, " beta_y: ", by
    print "q_x: ", qx, " q_y: ", qy
    print "beta_z: ", bz
    print "delta p/p: ", dpop

pr1 = Proton(energy)
pr1p = Proton(energy)
pr2 = Proton(energy)
pr2p = Proton(energy)
pr3 = Proton(energy)
pr3p = Proton(energy)
pr1.set_x(1.0e-4)
pr1p.set_npx(1.0e-4/bx)
pr2.set_y(1.0e-4)
pr2p.set_npy(1.0e-4/by)

bml.propagate(pr1)
bml.propagate(pr1p)
bml.propagate(pr2)
bml.propagate(pr2p)
if myrank == 0:
    print "propagated pr1: ", pr1.State()
    print "propagated pr1p: ", pr1p.State()
    print "propagated pr2: ", pr2.State()
    print "propagated pr2p: ", pr2p.State()

emitx /= (beta*gamma)
emity /= (beta*gamma)

if myrank == 0:
    print "Begin generating bunch..."

#bunch = synergia.optics.generate_matched_bunch_transverse(lattice_simulator, emitx, emity, stdz, dpop,
#                                        real_particles, macro_particles,
#                                        seed=seed)

bunch = synergia.optics.generate_matched_bunch(lattice_simulator, np.sqrt(emitx * bx), np.sqrt(emity * by), stdz, real_particles, macro_particles,
                                        seed=seed)

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
    print "expected stdz: ", stdz, " generated (on rank 0): ", np.std(particles[:,4])

if opts.spacecharge:
    space_charge = synergia.collective.Space_charge_3d_open_hockney(bunch.get_comm(), grid)
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, space_charge, num_steps)
else:
    no_op = synergia.simulation.Dummy_collective_operator("stub")
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator, no_op,num_steps)

multi_diagnostics_step = synergia.bunch.Multi_diagnostics()
# not doing step diagnostics
#multi_diagnostics_step.append(synergia.bunch.Diagnostics_full2(bunch, "mi_step_full2.h5"))

multi_diagnostics_turn = synergia.bunch.Multi_diagnostics()

if opts.turn_full2:
    multi_diagnostics_turn.append(synergia.bunch.Diagnostics_full2(bunch,"mi_full2.h5"))

if opts.turn_particles:
    multi_diagnostics_turn.append(synergia.bunch.Diagnostics_particles(bunch, "mi_particles.h5"))

# enable track saving
for part in range(0, opts.turn_tracks):
    multi_diagnostics_turn.append(synergia.bunch.Diagnostics_track(bunch,"turn_track_%02d.h5"%part,part))

if myrank == 0:
    print "Begin propagate..."

propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, num_turns, multi_diagnostics_step,
                     multi_diagnostics_turn, verbose)

if myrank == 0:
    print "Finish propagate"
