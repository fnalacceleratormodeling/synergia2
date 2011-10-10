#!/usr/bin/env python
import sys
import time
import synergia
import numpy as np
import re
import math
from synergia.optics.one_turn_map import linear_one_turn_map
from mu2e_options import opts
import mpi4py.MPI as MPI
from basic_toolkit import *
from mxyzptlk import *
from beamline import *
from physics_toolkit import *

###############################################################################
# Class for magnet ramping
###############################################################################
class Ramp_actions(synergia.simulation.Propagate_actions):
    def __init__(self, ramp_turns, turns_to_extract, final_k2l, initial_k1,
            final_k1):
        synergia.simulation.Propagate_actions.__init__(self)
        self.ramp_turns = ramp_turns
        self.turns_to_extract = turns_to_extract
        self.final_k2l = final_k2l
        self.initial_k1 = initial_k1
        self.final_k1 = final_k1
    def turn_end_action(self, stepper, bunch, turn_num):
        for element in stepper.get_lattice_simulator().get_lattice().get_elements():        
            # sextupole ramping...
            if turn_num <= ramp_turns:
                index = 0
                if element.get_type() == "multipole":
                    new_k2l = self.final_k2l[index] * turn_num / ramp_turns
                    element.set_double_attribute("k2l", new_k2l)
                    index += 1
                    if myrank == 0 and element.get_name() == "ddd_20_1" :
                        print
                        print "    turn                       :", turn_num
                        print "    updated thinsSextupole     :", element.get_name()
                        print "    initial k2l:               :", self.final_k2l[index - 1], "1/m^3"
                        print "    new k2l:                   :", new_k2l, "1/m^3"
            # quadrupole ramping...
            if turn_num > ramp_turns:
                index = 0
                epsilon = 1.0 * (turn_num - self.ramp_turns) / self.turns_to_extract
                if element.get_type() == "quadrupole":
                    new_k1 = (1.0 - epsilon) * self.initial_k1[index] + epsilon * self.final_k1[index]
                    element.set_double_attribute("k1", new_k1)
                    index += 1
                    if myrank == 0 and element.get_name() == "hqf1":
                        print
                        print "    turn                      :", turn_num
                        print "    updated quadrupole        :", element.get_name()
                        print "    epsilon                   :", epsilon
                        print "    initial k1                :", self.initial_k1[index - 1], "1/m^2"
                        print "    final k1                  :", self.final_k1[index - 1], "1/m^2"
                        print "    new k1                    :", new_k1, "1/m^2"
        stepper.get_lattice_simulator().update()

###############################################################################
#    adjust tune
###############################################################################
def adjust_tune(beamline_context, tune_h, tune_v, tolerance, tune_step):
    nu_h = beamline_context.getHorizontalFracTune()
    nu_v = beamline_context.getVerticalFracTune()
    final_tune_h = tune_h - int(tune_h)
    final_tune_v = tune_v - int(tune_v)

    if myrank == 0: 
        print "        Initial tunes: horizontal:", beamline_context.getHorizontalFracTune() 
        print "                       vertical  :", beamline_context.getVerticalFracTune() 
        print "        Final tunes:   horizontal:", final_tune_h
        print "                       vertical  :", final_tune_v

    dtune_h = (final_tune_h - nu_h) / (1 + int(np.abs(final_tune_h - nu_h)) / tune_step)
    target_tune_h = nu_h + dtune_h
    target_tune_v = final_tune_v

    step = 1
    while (np.abs(final_tune_h - nu_h) > 10.0 * tolerance):
        if myrank == 0:
            print "        Step", step
        while (np.abs(final_tune_h - nu_h) > tolerance):
            beamline_context.changeTunesBy(target_tune_h - nu_h, 0)
            #if ( errorCode != BeamlineContext::OKAY ) {
            #    if ( errorCode == BeamlineContext::NO_TUNE_ADJUSTER ) {
            #        (*outstreamptr) << "\n\n*** ERROR *** No tune adjuster in the context!" << endl;
            #    }
            #    exit(-1);
            #}
            beamline_context.changeTunesBy(0, target_tune_v - nu_v);
            #if ( errorCode != BeamlineContext::OKAY ) {
            #    if ( errorCode == BeamlineContext::NO_TUNE_ADJUSTER ) {
            #        (*outstreamptr) << "\n\n*** ERROR *** No tune adjuster in the context!" << endl;
            #    }
            #    exit(-1);
            #}
            nu_h = beamline_context.getHorizontalFracTune();
            nu_v = beamline_context.getVerticalFracTune();
            if myrank == 0:
                print "        Tunes: horizontal:", nu_h, "=", target_tune_h,
                print "+", nu_h - target_tune_h
                print "             : vertical  :", nu_v
        target_tune_h += dtune_h
        step += 1

###############################################################################
#   resonance sum
###############################################################################
def nint(x):
    lower = math.floor(x)
    upper = math.ceil(x)
    if  x - lower >= upper - x:
        return int(upper)
    else:
        return int(lower)

def resonance_sum(information, beamline, bhro):
    g = np.zeros([2],dtype=complex)
    imaginary_number = 1j
    factor = imaginary_number / (24.0 * np.sqrt(2.0) * np.pi * brho)

    num_elements = beamline.countHowMany()
    circumference = information[num_elements - 1].arcLength
    nu_x = information[num_elements - 1].psi.hor / (2.0 * np.pi)
    harmonic = float(nint(3.0 * nu_x))
    delta = nu_x - harmonic / 3.0

    if myrank == 0:
        print "        number of elements     :", num_elements
        print "        circumference          :", circumference, "m"
        print "        nu_x                   :", nu_x
        print "        harmonic               :", harmonic
        print "        delta                  :", delta

    index = 0
    for element in beamline:
        local_info = information[index]
        if element.Type() == "thinSextupole":
            theta = 2.0 * np.pi * local_info.arcLength / circumference
            beta_x = local_info.beta.hor
            psi_x = local_info.psi.hor
            phase = 3.0 * (psi_x - delta * theta)
            while np.pi < phase:
                phase -= 2.0 * np.pi
            while phase <= - np.pi:
                phase += 2.0 * np.pi
            if myrank == 0:
                print
                print "        index                  :", index
                print "        type                   :", element.Type()
                print "        name                   :", element.Name()
                print "        strength               :", element.Strength()
                print "        theta                  :", theta
                print "        beta_x                 :", beta_x
                print "        psi_x                  :", psi_x
                print "        phase                  :", phase
            beta32 = beta_x * np.sqrt(beta_x)
            increment = beta32 * 2.0 * element.Strength() * (np.cos(phase) - np.sin(phase) * imaginary_number)

            if re.search('ddd_50', element.Name()):
                g[0] += increment
            elif re.search('ddd_20', element.Name()):
                g[1] += increment
            else:
                print "Incorrectly named harmonic sextupole, perhaps?"
                sys.out(1)
        index += 1

    g[0] *= factor
    g[1] *= factor

    return g

###############################################################################
#   quick and dirty twiss parameter calculator from 2x2 Courant-Snyder map array
###############################################################################
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

###############################################################################
#   Mu2e Third Integer Resonant Extraction Simulation
###############################################################################

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
    print 
    print "Run Summary"
    print "    Macro particles                  :", macro_particles
    print "    Real particles                   :", real_particles
    print "    num_turns                        :", num_turns
    print "    num steps/turn                   :", num_steps
    print "    grid                             :", grid
    print "    map order                        :", map_order
    print "    (normalized) transverse emittance:", opts.norm_emit * 1e6, "mm-mrad"
    print "    stdz                             :", stdz
    print "    bunch length                     :", bunchlen_sec * 1e9, "nsec"
    print "    Radius aperture                  :", radius, "m"
    print "    RF Voltage                       :", rf_voltage, "MV"
    print "    X offset                         :", x_offset
    print "    Y offset                         :", y_offset
    print "    Z_offset                         :", z_offset
    print "    Space Charge is", 
    if opts.spacecharge:
        print "ON"
    else:
        print "not ON"
    print "    Generating diagnostics           :",
    if opts.turn_full2:
        print "turn_full2",
    if opts.turn_particles:
        print "turn_particles",
    if opts.turn_tracks:
        print "turn_tracks:", opts.turn_tracks,
    print

synergia_lattice = synergia.lattice.Mad8_reader().get_lattice("debunch", 
                "Debunch_modified_rf.lat")
synergia_elements = synergia_lattice.get_elements()

# with the aperture, all the particles are immediately eliminated
if radius > 0.0:
    for elem in synergia_elements:
        elem.set_double_attribute("aperture_radius", radius)

lattice_length = synergia_lattice.get_length()

reference_particle = synergia_lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
brho = reference_particle.get_momentum() / (synergia.foundation.pconstants.c / 1e9)
bunchlen_m = bunchlen_sec * beta * synergia.foundation.pconstants.c

emitx /= (beta * gamma) * 6.0
emity /= (beta * gamma) * 6.0

if myrank == 0:
    print "    lattice length                   :", lattice_length, "m"
    print "    brho                             :", brho, "T-m"
    print "    energy                           :", energy, "GeV"
    print "    beta                             :", beta
    print "    gamma                            :", gamma
    print "    bunch length (in second)         :", bunchlen_sec * 1e9, "nsec"
    print "    bunch length (in meter)          :", bunchlen_m, "m"

##############################################################################
#   rf cavity is not implemented yet (FIXME)
##############################################################################
# set rf cavity frequency (= harmno * beta * c / ring_length)
if myrank == 0:
    print
    print "Begin setting RF voltage..."
harmno = 4
freq = harmno * beta * synergia.foundation.pconstants.c/lattice_length
for element in synergia_elements:
    if element.get_type() == "rfcavity":
        element.set_double_attribute("volt", rf_voltage)
        element.set_double_attribute("freq", freq)
        if myrank == 0:
            print "    rfcavity                         :", element.get_name()
            print "    rf voltage                       :", element.get_double_attribute("volt"), "MV"
            print "    rf frequency                     :", element.get_double_attribute("freq") * 1e-6, "MHz"

###############################################################################
#   Set lattice_simulator and chef_lattice
###############################################################################
if myrank == 0:
    print
    print "Set lattice_simulator and chef_lattice"

lattice_simulator = synergia.simulation.Lattice_simulator(synergia_lattice, 
                map_order)
chef_lattice = lattice_simulator.get_chef_lattice()

###############################################################################
#   Set the lattice for third integer resonant extraction
###############################################################################
if myrank == 0:
    print
    print "Set the lattice for the third integer resonant extraction"
#   Establish the Jet environment and initiate a beamline context
    print
    print "....Establish the Jet environment and initiate a beamline context"
JetParticle.createStandardEnvironments(map_order)
probe = Proton(energy)
probe.setStateToZero()
beamline = chef_lattice.get_beamline()
beamline_context = BeamlineContext(probe, beamline) 

#   Set up the tune control circuit
if myrank == 0:
    print
    print "....Set the tune control circuit"
for element in beamline:
    if re.search('hqf1', element.Name()) or re.search('hqf2', element.Name()):
        beamline_context.addHTuneCorrectorQuadrupolePtr(element)
    elif re.search('hqd1', element.Name()) or re.search('hqd2', element.Name()):
        beamline_context.addVTuneCorrectorQuadrupolePtr(element)

#   Calculate target emittance from invariant emittance
#   (This is the area of a triangle tangent to a circle.)
if myrank == 0:
    print
    print "....Calculate target emittance from invariant emittance"
safety_factor = 1.05
emittance_x = opts.norm_emit * np.pi / (beta * gamma)
emittance_x *= safety_factor 
emittance_x *= 3.0 * np.sqrt(3.0) / np.pi

#   Set initial tunes.
if myrank == 0:
    print
    print "....Set initial tunes"
bare_tune = [opts.tuneh, opts.tunev]

#   Adjust Tunes
if myrank == 0:
    print
    print "....Adjust Tunes"
tune_tolerance = 1.0e-8
adjuster_tune_step = 0.005
adjust_tune(beamline_context, bare_tune[0], bare_tune[1], tune_tolerance, adjuster_tune_step)

#   Store initial quad settings for initial tune
if myrank == 0:
    print
    print "....Store initial quad settings for initial tune"
initial_k1 = []
index = 0
for element in beamline:
    if element.Type() == "quadrupole":
        initial_k1.append(element.Strength() / brho)       # 1/m^2
        index += 1
        #if myrank == 0:
        #    print "       ", element.Name(), initial_k1[index - 1], "=", element.Strength(), "/ brho"

#   Determine the target value of |g|
if myrank == 0:
    print
    print "....Determine the target value of |g|"
resonant_tune = opts.resonant_tune
delta = bare_tune[0] - resonant_tune
target_abs_g = np.abs(delta) / np.sqrt(2.0 * np.sqrt(3.0) * emittance_x)
target_phase_g = opts.phase_g * np.pi / 180.0
target_g = target_abs_g * (np.cos(target_phase_g) + 1j * np.sin(target_phase_g))
if myrank == 0:
    print "        resonant tune          :", resonant_tune
    print "        target |g|             :", target_abs_g
    print "        target psi_g           :", opts.phase_g, "degree"

#   Calculate current value of |g|
if myrank == 0:
    print
    print "....Calculate current value of g"
    print "........Calling resonance_sum"
    print 
information = beamline_context.getTwissArray()
g = np.zeros([2], dtype=complex)
g = resonance_sum(information, beamline, brho)
total_g = g[0] + g[1]
abs_g = np.abs(total_g)

if myrank == 0:
    print 
    print "        current value: g       :", total_g
    print "                      |g|      :", abs_g

#   Adjust harmonic sextupole strengths
if myrank == 0:
    print
    print "....Adjust harmonic sextupole strengths: Setting (integrated) sextupole strengths"
    print "........Change polarity, if desired"
if opts.flip:
    for element in beamline:
        if element.Type() == "thinSextupole":
            element.setStrength(- element.Strength())

ratio = target_abs_g / abs_g
g *= ratio
if myrank == 0:
    print
    print "........Scale to target value of harmonic coupling: Setting (integrated) sextupole strengths"
    print "            ratio              :", ratio
for element in beamline:
    if element.Type() == "thinSextupole":
        old_strength = element.Strength()
        new_strength = element.Strength() * ratio
        element.setStrength(new_strength)
        if myrank == 0:
            print 
            print "            type               :", element.Type()
            print "            name               :", element.Name()
            print "            strength           :", old_strength, "->", new_strength, "T/m"
            print "                               : B''l =", 2.0 * new_strength, "T/m"

if myrank == 0:
    print 
    print "........Rotate, if desired"
denom = (np.conjugate(g[1]) * g[0]).imag
s = np.zeros([2], 'd')
s[0] = (np.conjugate(g[1]) * target_g).imag / denom;
s[1] = -(np.conjugate(g[0]) * target_g).imag / denom;
if myrank == 0:
    print "            g                  :", target_g
    print "            g1                 :", g[0]
    print "            g2                 :", g[1]
    print "            s1                 :", s[0]
    print "            s2                 :", s[1]
    print "            target_g           :", np.abs(target_g)
    print "            g                  :", np.abs(g[0] + g[1])
    print "            target_phase_g     :", target_phase_g
    print "            delta              :", delta
    print
    print "............Doing rotation"
g[0] *= s[0]
g[1] *= s[1]
for element in beamline:
    if element.Type() == "thinSextupole":
        if re.search('ddd_50', element.Name()):
            old_strength = element.Strength()
            new_strength = old_strength * s[0]
            element.setStrength(new_strength)
        elif re.search('ddd_20', element.Name()):
            old_strength = element.Strength()
            new_strength = old_strength * s[1]
            element.setStrength(new_strength)
        if myrank == 0:
            print
            print "            type               :", element.Type()
            print "            name               :", element.Name()
            print "            strength           :", old_strength, "->", new_strength, "T/m"
            print "                               : B''l =", 2.0 * new_strength, "T/m"

#   Store final settings of harmonic sextupole circuits
final_k2l = []
index = 0
if myrank == 0:
    print
    print "........Store final settings of harmonic sextupole circuits"
for element in beamline:
    if element.Type() == "thinSextupole":
        final_k2l.append(element.Strength() / brho)        # 1/m^3
        index += 1
        if myrank == 0:
            print 
            print "            type               :", element.Type()
            print "            name               :", element.Name()
            print "            strength (k2l)     :", final_k2l[index-1], "1/m^3"

#   Set resonant tune
if myrank == 0:
    print
    print "....Set resonant tune"
tune_h = resonant_tune
adjust_tune(beamline_context, tune_h, bare_tune[1], tune_tolerance, adjuster_tune_step)

#   Store final quad settings for resonant tune
if myrank == 0:
    print
    print "....Store final quad settings for resonat tune"
final_k1 = []
index = 0
for element in beamline:
    if element.Type() == "quadrupole":
        final_k1.append(element.Strength() / brho)         # 1/m^2
        index += 1
        #if myrank == 0:
        #    print "       ", element.Name(), final_k1[index - 1], "=", element.Strength(), "/ brho"

#   Testing
if myrank == 0:
    print
    print "....Test: after setting final conditions"
index = 0
for element in beamline:
    if element.Type() == "quadrupole":
        element.setStrength(final_k1[index])
        index += 1
#beamline_context.reset()
nu_h = beamline_context.getHorizontalEigenTune()
nu_v = beamline_context.getVerticalEigenTune()
fractional_tune = resonant_tune - int(resonant_tune)
if myrank == 0:
    print
    print "        Fractional tunes: horizontal:", nu_h, "= resFracTune +", 
    print nu_h - fractional_tune
    print "                          vertical  :", nu_v
    print
    print "....Test: after setting initial conditions"
index = 0
for element in beamline:
    if element.Type() == "quadrupole":
        element.setStrength(initial_k1[index])
        index += 1
#beamline_context.reset()
nu_h = beamline_context.getHorizontalEigenTune()
nu_v = beamline_context.getVerticalEigenTune()
fractional_tune = bare_tune[0] - int(bare_tune[0])
resFracTune = resonant_tune - int(resonant_tune)
if myrank == 0:
    print
    print "        Fractional tunes: horizontal:", nu_h, "= resFracTune +",
    print nu_h - fractional_tune
    print "                          vertical  :", nu_v

#   Determine ramping turns
ramp_turns = opts.rampturns
extraction_fraction = opts.extraction_fraction
turns_to_extract = int(extraction_fraction * beta * synergia.foundation.pconstants.c / lattice_length / 60.0 ) + ramp_turns

if myrank == 0:
    print
    print "....Determine ramping turns"
    print "        sextupole ramping turns:", ramp_turns
    print "        extraction fraction    :", extraction_fraction
    print "        turns to extract:      :", turns_to_extract
    print "        final turn             :", num_turns

ramp_actions = Ramp_actions(ramp_turns, turns_to_extract, final_k2l,
                initial_k1, final_k1)

###############################################################################
#   Set initial element strengths for beam matching
###############################################################################
if myrank == 0:
    print
    print "Set initial element strengths"
    print
    print "....Set initial quad strengths of the Synergia lattice"
index = 0
for element in synergia_elements:
    if element.get_type() == "quadrupole":
        element.set_double_attribute("k1", initial_k1[index])
        index += 1
if myrank == 0:
    print
    print "....Set initial quad strengths of the CHEF lattice"
index = 0
for element in beamline:
    if element.Type() == "quadrupole":
        element.setStrength(initial_k1[index] * brho)
        index += 1

if myrank == 0:
    print
    print "....Set initial sextupole strengths of the CHEF lattice"
index = 0
for element in beamline:
    if element.Type() == "thinSextupole":
        element.setStrength(0.0)
        index += 1

###############################################################################
#  CHEF One Turn Map 
###############################################################################
if myrank == 0:
    print
    print "CHEF One Turn Map"
JetParticle.createStandardEnvironments(map_order)
jpr = JetProton(energy)
#bml = chef_lattice.get_beamline()
#bml.propagate(jpr)
beamline.propagate(jpr)
chefmap = jpr.State()
chefoneturnmap = chefmap.jacobian()
if myrank == 0:
    print chefoneturnmap

###############################################################################
#   One Turn Map from Synergia
###############################################################################
map = linear_one_turn_map(lattice_simulator)
if myrank == 0:
    print "One turn map from synergia2.5 infrastructure"
    print np.array2string(map, max_line_width=200, precision=3)
    print "    det(M) :", np.linalg.det(map)

# checking eigen vector and values of map
#[l, v] = np.linalg.eig(map)
#print "l:", l
#print "v:", v
#if myrank == 0:
#    print "    eigenvalues:"
#    for z in l: 
#        print "    |z|:", abs(z), " z:", z, " tune:", np.log(z).imag/(2.0*np.pi)

[ax, bx, qx] = map2twiss(map[0:2, 0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6, 4:6])
alpha_twiss, beta_twiss = synergia.optics.get_alpha_beta(map)
if ax != alpha_twiss[0] or ay != alpha_twiss[1] or bx != beta_twiss[0] or by != beta_twiss[1]:
    print "Error:"    #FIXME
    sys.out(1)
dpop = bunchlen_m/bz

if myrank == 0:
    print
    print "Lattice parameters (assuming uncoupled map)"
    print "    alpha_x                          :", alpha_twiss[0]
    print "    alpha_y                          :", alpha_twiss[1]
    print "    beat_x                           :", beta_twiss[0]
    print "    beta_y                           :", beta_twiss[1]
    print "    q_x                              :", qx
    print "    q_y                              :", qy
    print "    beta_z                           :", bz
    print "    delta_p/p                        :", dpop

emittance = [emitx, emity]
width, r = synergia.optics.match_transverse_twiss_emittance(emittance,
                alpha_twiss, beta_twiss)
if myrank == 0:
    print
    print "Beam parameters for matching"
    print "    emitx                            :", emitx * 1e6, "mm-mrad"
    print "    emity                            :", emity * 1e6, "mm-mrad"
    print "    xwidth                           :", np.sqrt(emitx * bx), "m"
    print "    ywidth                           :", np.sqrt(emity * by), "m"

if myrank == 0:
    print
    print "Begin generating bunch..."
#bunch = synergia.optics.generate_matched_bunch_transverse(lattice_simulator, 
#                emitx, emity, stdz, dpop, real_particles, macro_particles,
#                seed=seed)
bunch = synergia.optics.generate_matched_bunch(lattice_simulator,
                np.sqrt(emitx * beta_twiss[0]), np.sqrt(emity * beta_twiss[1]),
                bunchlen_m, real_particles, macro_particles, seed=seed)

# apply offset to bunch
particles = bunch.get_local_particles()

particles[:,0] = particles[:,0] + x_offset
particles[:,2] = particles[:,2] + y_offset
particles[:,4] = particles[:,4] + z_offset

if myrank == 0:
    print
    print "    expected stdx:", np.sqrt(emitx * bx), "(m) generated (on rank 0):", np.std(particles[:,0]), "(m)"
    print "    expected stdy:", np.sqrt(emity * by), "(m) generated (on rank 0):", np.std(particles[:,2]), "(m)"
    print "    expected stdz:", bunchlen_m, "   (m) generated (on rank 0):", np.std(particles[:,4]), "   (m)"

###############################################################################
#   Collective operator
###############################################################################
solver = opts.spacecharge
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
diagnostics_actions = synergia.simulation.Standard_diagnostics_actions()
#diagnostics per step
for part in range(0, opts.step_tracks):
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_track(bunch,
                            "mu2e_step_track_%02d.h5" % part, part))
if opts.step_full2:
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_full2(bunch,
                            "mu2e_step_full2.h5"))
if opts.step_particles:
    diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_particles(bunch,
                            "mu2e_step_particles.h5"))
# diagnostics per turn
for part in range(0, opts.turn_tracks):
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_track(bunch,
                            "mu2e_turn_track_%02d.h5" % part, part))
if opts.turn_full2:
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_full2(bunch,
                            "mu2e_turn_full2.h5"))
if opts.turn_particles:
    diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_particles(bunch,
                            "mu2e_turn_particles.h5"))

###############################################################################
#   Propagate and track particles
###############################################################################
if myrank == 0:
    print
    print "Begin propagate..."

t0 = time.time()
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, num_turns, diagnostics_actions, ramp_actions, 
                verbose)
t1 = time.time()

if myrank == 0:
    print
    print "propagate time =", t1 - t0
