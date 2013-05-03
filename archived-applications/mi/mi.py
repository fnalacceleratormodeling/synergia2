#!/usr/bin/env python
import sys
import synergia
import numpy as np
from mi_options import opts
import mpi4py.MPI as MPI

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
    print "(normalized) transverse emittance: ", opts.norm_emit
    print "stdz: ", stdz
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

lattice = synergia.lattice.Mad8_reader().get_lattice("ring_p_q605", "mi20-egs-thinrf.lat")
lattice_length = lattice.get_length()

# set all elements to use chef_propagate
# and MI elliptical aperture
if myrank == 0:
    print "Begin constructing lattice... "
for elem in lattice.get_elements():
    elem.set_string_attribute("extractor_type", "chef_propagate")

    # these aperture sizes come from _Main_Injector_Transverse_Apertures_,
    # Bruce C. Brown, 03/13/2006.
    elem.set_string_attribute("aperture_type", "elliptical")
    elem.set_double_attribute("elliptical_aperture_horizontal_radius", 60.26/1000.0)
    elem.set_double_attribute("elliptical_aperture_vertical_radius", 23.69/1000.0)

if myrank == 0:
    print "Finished constructing lattice"


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

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)

map = lattice_simulator.get_linear_one_turn_map()
if myrank==0:
    print "one turn map from synergia2.5 infrastructure"
    print np.array2string(map, max_line_width=200)


[l, v] = np.linalg.eig(map)

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

emitx /= (beta*gamma)
emity /= (beta*gamma)

if myrank == 0:
    print "Begin generating bunch..."

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


bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

if opts.turn_full2:
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("mi_full2.h5"))

if opts.turn_particles:
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("mi_particles.h5"))

# enable track saving for the first 100 turns
if  opts.turn_tracks:
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("turn_track.h5", opts.turn_tracks), range(100))

if myrank == 0:
    print "Begin propagate..."

propagator = synergia.simulation.Propagator(stepper)
max_turns = 0
propagator.propagate(bunch_simulator, num_turns, max_turns, verbose)

if myrank == 0:
    print "Finish propagate"
