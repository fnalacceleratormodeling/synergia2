#!/usr/bin/env python
import sys
import os
import numpy as np
import synergia

from foborodobo128_options import opts

#####################################

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        raise RuntimeError, "map is unstable"

    mu =np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

#######################################################
################################################################################

def print_bunch_stats(bunch, fo):
    coord_names = ("x", "xp", "y", "yp", "c*dt", "dp/p")

    means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
    stds = synergia.bunch.Core_diagnostics().calculate_std(bunch, means)
    print >>fo, "%20s   %20s   %20s"%("coord","mean","rms")
    print >>fo, "%20s   %20s   %20s"%("====================",
                                      "====================",
                                      "====================")
    for i in range(6):
        print >>fo, "%20s   %20.12e   %20.12e"%(coord_names[i], means[i], stds[i]\
)


################################################################################

logger = synergia.utils.Logger(0)

lattice = synergia.lattice.MadX_reader().get_lattice("model", "foborodobo128.madx")

refpart = lattice.get_reference_particle()
energy = refpart.get_total_energy()
momentum = refpart.get_momentum()
gamma = refpart.get_gamma()
beta = refpart.get_beta()

print >>logger, "energy: ", energy
print >>logger, "momentum: ", momentum
print >>logger, "gamma: ", gamma
print >>logger, "beta: ", beta

length = lattice.get_length()
print >>logger, "lattice length: ", length

harmon = 32.0
freq = harmon * beta * synergia.foundation.pconstants.c/length
print >>logger, "frequency: ", freq

for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("freq", freq*1.0e-6)
        print >>logger, "rfcavity: ", elem.print_()

for elem in lattice.get_elements():
    elem.set_string_attribute("extractor_type", "chef_propagate")

for elem in lattice.get_elements():
    elem.set_string_attribute("extractor_type", "libff")

f = open("foborodobo128_lattice.out", "w")
print >>f, lattice.as_string()
f.close()

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)

myrank = 0
map = lattice_simulator.get_linear_one_turn_map()
print >>logger, "one turn map from synergia2.5 infrastructure"
print >>logger, np.array2string(map, max_line_width=200)


[l, v] = np.linalg.eig(map)

#print "l: ", l
#print "v: ", v

print >>logger, "eigenvalues: "
for z in l:
    print >>logger, "|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi)

[ax, bx, qx] = map2twiss(map[0:2,0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6,4:6])

print >>logger, "Lattice parameters (assuming uncoupled map)"
print >>logger, "alpha_x: ", ax, " alpha_y: ", ay
print >>logger, "beta_x: ", bx, " beta_y: ", by
print >>logger, "q_x: ", qx, " q_y: ", qy
print >>logger, "q_z: ", qz, " beta_z: ", bz


#lattice_simulator.print_lattice_functions()

alpha_c = lattice_simulator.get_momentum_compaction()
slip_factor = alpha_c - 1/gamma**2
print >>logger, "alpha_c: ", alpha_c, ", slip_factor: ", slip_factor

hchrom = lattice_simulator.get_horizontal_chromaticity()
vchrom = lattice_simulator.get_vertical_chromaticity()

print >>logger, "horizontal chromaticity: %.16g"%hchrom
print >>logger, "vertical chromaticity: %.16g"%vchrom

chef_beamline = lattice_simulator.get_chef_lattice().get_beamline()
f = open("foborodobo128_beamline.out","w")
print >>f, synergia.lattice.chef_beamline_as_string(chef_beamline)
f.close()

stdx = opts.stdx
stdy = opts.stdy
stdz = opts.stdz

print >>logger, "stdx: ", stdx
print >>logger, "stdy: ", stdy
print >>logger, "stdz: ", stdz

macro_particles = opts.macro_particles
real_particles = opts.real_particles
print >>logger, "macro_particles: ", macro_particles
print >>logger, "real_particles: ", real_particles

if opts.matching == "6dmoments":
    print >>logger, "Matching with 6d moments"
    bunch = synergia.optics.generate_matched_bunch(lattice_simulator, stdx, stdy, stdz, real_particles, macro_particles, rms_index=[0,2,4], periodic=False)
else:
    print >>logger, "Matching with normal form"
    actions = lattice_simulator.get_stationary_actions(stdx, stdy, stdz/beta)
    comm = synergia.utils.Commxx()
    bunch = synergia.bunch.Bunch(refpart, macro_particles, real_particles, comm)
    seed = opts.seed
    dist = synergia.foundation.Random_distribution(seed, comm)
    synergia.simulation.populate_6d_stationary_gaussian(dist, bunch, actions, lattice_simulator)

print_bunch_stats(bunch, logger)

# laslett tune shift
#
# Delta_Q = N * r_0 * 2.0 * F*G/ (pi * beta**2 * gamma**3 * emittance95 *B)
#
# F=laslett coefficient=0.5 for circular beam
# G=Form factor for transverse distribution = 2.0 for gaussian beam
# B=Bunching factor = mean/peak longitudinal density ~0.3
# emittance95 = 4*sigma**2/beta
#

emit_x = 4.0*stdx**2/bx
F=0.5
G=2.0
Bf=0.3
laslett = opts.real_particles*F*G*synergia.foundation.pconstants.rp/(np.pi*beta**2*gamma**3*emit_x*Bf)

print >>logger, "Laslett tune shift: ", laslett

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("full2.h5"))
if opts.particles:
    print >>logger, "saving particles"
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("particles.h5"))
if opts.tracks:
    print >>logger, "saving ", opts.tracks, " particle tracks"
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("tracks.h5", opts.tracks))

steppertype = opts.stepper
if opts.spacecharge and steppertype != "splitoperator":
    print >>logger, "changing stepper to splitoperator because spacecharge is ON"
    steppertype = "splitoperator"

if steppertype == "independent":
    print >>logger, "using independent stepper ", opts.steps, " steps/turn"
elif steppertype == "elements":
    print >>logger, "using steps by element ", opts.steps, " steps/element"
elif steppertype == "splitoperator":
    print >>logger, "using split operator stepper ", opts.steps, " steps/turn"
else:
    raise  RuntimeError, "unknown stepper type: %s"%steppertype

print "rank ", comm.get_rank(), " propagating bunch with ", bunch.get_local_num(), " particles"

if opts.spacecharge:
    print >>logger, "space charge is ON"
    # space charge requires a communicator for collective effects
    coll_comm = synergia.utils.Commxx(True)
    print >>logger, "space charge grid: ", opts.gridx, opts.gridy, opts.gridz
    grid = [opts.gridx, opts.gridy, opts.gridz]
    coll_operator = synergia.collective.Space_charge_2d_open_hockney(coll_comm, grid)
else:
    coll_operator = synergia.simulation.Dummy_collective_operator("stub")
    print >>logger, "space charge is OFF"

if steppertype == "independent":
    stepper = synergia.simulation.Independent_stepper(lattice, 1, opts.steps)
elif steppertype == "elements":
    stepper = synergia.simulation.Independent_stepper_elements(lattice, 1, opts.steps)
elif steppertype == "splitoperator":
    stepper = synergia.simulation.Split_operator_stepper(lattice, 1, coll_operator, opts.steps)

propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(100)

#propagator.propagate(bunch_simulator, 100, 100, 1)
# propagator.propagate(bunch_simulator, turns, max_turns, verbosity)
propagator.propagate(bunch_simulator, opts.turns, opts.turns, 1)
