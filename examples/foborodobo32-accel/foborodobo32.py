#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import numpy as np
import synergia

from foborodobo32_options import opts
from ramp_module import Ramp_actions

#####################################

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        raise RuntimeError("map is unstable")

    mu =np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

#######################################################
# get focussing and defocussing quadrupoles for adjust tune
def get_fd_quads(lattice):
    f_quads = []
    d_quads = []
    for elem in lattice.get_elements():
        if elem.get_type() == "quadrupole":
            if elem.get_double_attribute("k1") > 0.0:
                f_quads.append(elem)
            elif elem.get_double_attribute("k1") < 0.0:
                d_quads.append(elem)
    return (f_quads, d_quads)

################################################################################

def print_bunch_stats(bunch, fo):
    coord_names = ("x", "xp", "y", "yp", "c*dt", "dp/p")

    means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
    stds = synergia.bunch.Core_diagnostics().calculate_std(bunch, means)
    print("%20s   %20s   %20s"%("coord","mean","rms"), file=fo)
    print("%20s   %20s   %20s"%("====================",
                                      "====================",
                                      "===================="), file=fo)
    for i in range(6):
        print("%20s   %20.12e   %20.12e"%(coord_names[i], means[i], stds[i]), file=fo)


################################################################################

numbers = [3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3]

logger = synergia.utils.Logger(0)

# read the lattice in from a MadX sequence file
lattice = synergia.lattice.MadX_reader().get_lattice("model", "foborodobo32.madx")

# set the momentum of the beam.  This could have been in a beam statement
# in the lattice file, but this gives the option of setting it on the
# command line by creating a reference particle with the desired momentum.

# create the Reference particle object with the correct momentum
momentum = opts.momentum
four_momentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.mp)
four_momentum.set_momentum(momentum)

refpart = synergia.foundation.Reference_particle(1, four_momentum)

# set it into the lattice object
lattice.set_reference_particle(refpart)

energy = refpart.get_total_energy()
momentum = refpart.get_momentum()
gamma = refpart.get_gamma()
beta = refpart.get_beta()

print("energy: ", energy, file=logger)
print("momentum: ", momentum, file=logger)
print("gamma: ", gamma, file=logger)
print("beta: ", beta, file=logger)

length = lattice.get_length()
print("lattice length: ", length, file=logger)

# set RF parameters. The RF voltage and phase (lag) are set to give a
# synchrotron tune and a stable bucket.
harmon = 32.0
freq = harmon * beta * synergia.foundation.pconstants.c/length

rf_volt = opts.rf_volt
rf_lag = opts.lag

print("RF frequency: ", freq, " Hz", file=logger)
print("RF voltage: ", rf_volt, " MV", file=logger)

for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", rf_volt)
        elem.set_double_attribute("lag", rf_lag)
        elem.set_double_attribute("freq", freq*1.0e-6)
        print("rfcavity: ", elem.print_(), file=logger)

f = open("foborodobo32_lattice.out", "w")
print(lattice.as_string(), file=f)
f.close()

# the lattice_simulator object lets us do some computations for
# lattice functions and other parameters.
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)

myrank = 0
map = lattice_simulator.get_linear_one_turn_map()
print("one turn map from synergia2.5 infrastructure", file=logger)
print(np.array2string(map, max_line_width=200), file=logger)

[l, v] = np.linalg.eig(map)

#print( "l: ", l)
#print( "v: ", v)

print("eigenvalues: ", file=logger)
for z in l:
    print("|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi), file=logger)

[ax, bx, qx] = map2twiss(map[0:2,0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6,4:6])

print("Lattice parameters (assuming uncoupled map)", file=logger)
print("alpha_x: ", ax, " alpha_y: ", ay, file=logger)
print("beta_x: ", bx, " beta_y: ", by, file=logger)
print("q_x: ", qx, " q_y: ", qy, file=logger)
print("q_z: ", qz, " beta_z: ", bz, file=logger)


#lattice_simulator.print_lattice_functions()

alpha_c = lattice_simulator.get_momentum_compaction()
slip_factor = alpha_c - 1/gamma**2
print("alpha_c: ", alpha_c, ", slip_factor: ", slip_factor, file=logger)

f_quads, d_quads = get_fd_quads(lattice)
print("len(f_quads): ", len(f_quads), " len(d_quads): ", len(d_quads), file=logger)

(orig_xtune, orig_ytune) = lattice_simulator.get_both_tunes()
print("Original base tunes, x: ", orig_xtune, " y: ", orig_ytune, file=logger)

# adjust tunes if requested by the xtune and/or the ytune parameter using
# the list of focussing or defocussing quadruples as adjusters.

do_adjust_tunes = False
if opts.xtune or opts.ytune:
    do_adjust_tunes = True
    if opts.xtune:
        target_xtune = opts.xtune
    else:
        target_xtune = orig_xtune
    if opts.ytune:
        target_ytune = opts.ytune
    else:
        target_ytune = orig_ytune

if do_adjust_tunes:
    print("adjusting tunes, x: ", opts.xtune," y: ", opts.ytune, file=logger)
    lattice_simulator.adjust_tunes(target_xtune, target_ytune, f_quads, d_quads, 1.0e-6, 1)

hchrom = lattice_simulator.get_horizontal_chromaticity()
vchrom = lattice_simulator.get_vertical_chromaticity()

print("horizontal chromaticity: %.16g"%hchrom, file=logger)
print("vertical chromaticity: %.16g"%vchrom, file=logger)

chef_beamline = lattice_simulator.get_chef_lattice().get_beamline()
f = open("foborodobo32_beamline.out","w")
print(synergia.lattice.chef_beamline_as_string(chef_beamline), file=f)
f.close()

stdx = opts.stdx
stdy = opts.stdy
stdz = opts.stdz

print("stdx: ", stdx, file=logger)
print("stdy: ", stdy, file=logger)
print("stdz: ", stdz, file=logger)

macro_particles = opts.macro_particles
real_particles = opts.real_particles
print("macro_particles: ", macro_particles, file=logger)
print("real_particles: ", real_particles, file=logger)

comm = synergia.utils.Commxx()

# generate a 6D matched bunch using either normal forms or a 6D moments procedure
if opts.matching == "6dmoments":
    print("Matching with 6d moments", file=logger)
    bunch = synergia.optics.generate_matched_bunch(lattice_simulator, stdx, stdy, stdz, real_particles, macro_particles, rms_index=[0,2,4], periodic=False)
else:
    print("Matching with normal form", file=logger)
    actions = lattice_simulator.get_stationary_actions(stdx, stdy, stdz/beta)
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

print("Laslett tune shift: ", laslett, file=logger)

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

# define the bunch diagnostics to save

bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("full2.h5"))
if opts.particles:
    print("saving particles", file=logger)
    # save_particles option=n, 0: don't save, non-zero: save n particles
    # particles_period=n option, save particles every n turns
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("particles.h5", 0, opts.particles), opts.particles_period)
if opts.tracks:
    print("saving ", opts.tracks, " particle tracks", file=logger)
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("tracks.h5", opts.tracks))


# define the stepper, propagator and run the simulation.

steppertype = opts.stepper
if opts.spacecharge and steppertype != "splitoperator":
    print("changing stepper to splitoperator because spacecharge is ON", file=logger)
    steppertype = "splitoperator"

if opts.spacecharge is None:
    spacecharge = None
    print("space charge is OFF", file=logger)
elif ((opts.spacecharge == "off") or (opts.spacecharge == "0")):
    spacecharge = None
    print("space charge is OFF", file=logger)
elif opts.spacecharge == "2d-bassetti-erskine":
    print("space charge 2d-bassetti-erskine is ON", file=logger)
    coll_operator = synergia.collective.Space_charge_2d_bassetti_erskine()
    coll_operator.set_longitudinal(0)
elif opts.spacecharge == "2d-openhockney":
    print("space charge 2d-openhockney is ON", file=logger)
    # openhockney space charge requires a communicator for collective effects
    coll_comm = synergia.utils.Commxx(True)
    print("space charge grid: ", opts.gridx, opts.gridy, opts.gridz, file=logger)
    grid = [opts.gridx, opts.gridy, opts.gridz]
    coll_operator = synergia.collective.Space_charge_2d_open_hockney(coll_comm, grid)
elif opts.spacecharge == "3d-openhockney":
    print("space charge 3d-openhockney is ON", file=logger)
    coll_comm = synergia.utils.Commxx(True)
    print("space charge grid: ", opts.gridx, opts.gridy, opts.gridz, file=logger)
    grid = [opts.gridx, opts.gridy, opts.gridz]
    coll_operator = synergia.collective.Space_charge_3d_open_hockney(coll_comm, grid)
else:
    raise RuntimeError("unknown space charge operator")


if steppertype == "independent":
    print("using independent stepper ", opts.steps, " steps/turn", file=logger)
elif steppertype == "elements":
    print("using steps by element ", opts.steps, " steps/element", file=logger)
elif steppertype == "splitoperator":
    print("using split operator stepper ", opts.steps, " steps/turn", file=logger)
else:
    raise  RuntimeError("unknown stepper type: %s"%steppertype)

print(comm.get_rank(), " propagating bunch with ", bunch.get_local_num(), " particles", file=logger)


if steppertype == "independent":
    stepper = synergia.simulation.Independent_stepper(lattice, 1, opts.steps)
elif steppertype == "elements":
    stepper = synergia.simulation.Independent_stepper_elements(lattice, 1, opts.steps)
elif steppertype == "splitoperator":
    stepper = synergia.simulation.Split_operator_stepper(lattice, 1, coll_operator, opts.steps)

ramp_actions = Ramp_actions(numbers)

propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(100)

#propagator.propagate(bunch_simulator, 100, 100, 1)
# propagator.propagate(bunch_simulator, turns, max_turns, verbosity)
propagator.propagate(bunch_simulator, ramp_actions, opts.turns, opts.turns, 1)
