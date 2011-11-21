#!/usr/bin/env python
import sys
import time
import synergia
import numpy
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
#   resonance sum
###############################################################################
def nint(x):
    lower = math.floor(x)
    upper = math.ceil(x)
    if  x - lower >= upper - x:
        return int(upper)
    else:
        return int(lower)

def resonance_sum(lattice_simulator, brho):
    imaginary_number = 1j
    factor = imaginary_number / ((6.0 * numpy.sqrt(2.0)) * (4.0 * numpy.pi) * brho)

    synergia_elements = lattice_simulator.get_lattice().get_elements()
    num_elements = len(synergia_elements)
    last_element = synergia_elements[num_elements-1]
    last_lattice_functions = lattice_simulator.get_lattice_functions(last_element)
    nu_x = last_lattice_functions.psi_x / (2.0 * numpy.pi)
    circumference = last_lattice_functions.arc_length
    harmonic = float(nint(3.0 * nu_x))
    delta = nu_x - harmonic / 3.0

    if myrank == 0:
        print "        number of elements           :", num_elements
        print "        nu_x                         :", nu_x
        print "        circumference                :", circumference, "m"
        print "        harmonic                     :", harmonic
        print "        delta                        :", delta

    g = numpy.zeros([2],dtype=complex)
    index = 0
    for element in synergia_elements:
#lattice_simulator.get_lattice().get_elements():
        if element.get_type() == "multipole":
            lattice_functions = lattice_simulator.get_lattice_functions(element)
            theta = 2.0 * numpy.pi * lattice_functions.arc_length / circumference
            beta_x = lattice_functions.beta_x
            beta32 = beta_x * numpy.sqrt(beta_x)
            psi_x = lattice_functions.psi_x
            phase = 3.0 * (psi_x - delta * theta)
            while numpy.pi < phase:
                phase -= 2.0 * numpy.pi
            while phase <= - numpy.pi:
                phase += 2.0 * numpy.pi
            if myrank == 0:
                print
                print "        index                        :", index
                print "        type                         :", 
                print element.get_type()
                print "        name                         :", 
                print element.get_name()
                print "        strength                     :", 
                print element.get_double_attribute("k2l")
                print "        strength in chef unit        :",
                print element.get_double_attribute("k2l") * brho / 2.0
                print "        theta                        :", theta
                print "        beta_x                       :", beta_x
                print "        psi_x                        :", psi_x
                print "        phase                        :", phase
            increment = beta32 * brho \
                    * element.get_double_attribute("k2l") \
                    * (numpy.cos(phase) - numpy.sin(phase) * imaginary_number)

            if re.search('ddd_50', element.get_name()):
                g[0] += increment
            elif re.search('ddd_20', element.get_name()):
                g[1] += increment
            else:
                print "Incorrectly named harmonic sextupole, perhaps?"
                sys.exit(1)
        index += 1
    g[0] *= factor
    g[1] *= factor
    return g
###############################################################################
#   quick and dirty twiss parameter calculator from 2x2 Courant-Snyder map 
#   array by Eric
###############################################################################
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0] + csmap[1,1])
    asinmu = 0.5 * (csmap[0,0] - csmap[1,1])

    if abs(cosmu) > 1.0:
        if myrank == 0:
            print "error, map is unstable"
    mu = numpy.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * numpy.pi - mu

    beta = csmap[0,1]/numpy.sin(mu)
    alpha = asinmu/numpy.sin(mu)
    tune = mu/(2.0*numpy.pi) 

    return (alpha, beta, tune)

###############################################################################
###############################################################################
#                                                                             #
#   Mu2e Third Integer Resonant Extraction Simulation                         #
#                                                                             #
###############################################################################
###############################################################################

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
emitx = opts.norm_emit
emity = opts.norm_emit
stdz = opts.stdz
bunchlen_sec = opts.bunchlen * 1e-9
rf_voltage = opts.rf_voltage

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
    print "    (normalized) transverse emittance:", opts.norm_emit * 1e6, 
    print "mm-mrad"
    print "    stdz                             :", stdz
    print "    bunch length                     :", bunchlen_sec * 1e9, "nsec"
    print "    Radius aperture                  :", radius, "m"
    print "    RF Voltage                       :", rf_voltage, "MV"
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

synergia_lattice = synergia.lattice.Mad8_reader().get_lattice("debunch",
                "Debunch_modified_rf.lat")
synergia_elements = synergia_lattice.get_elements()

# with the aperture, all the particles are immediately eliminated
if radius > 0.0:
    for element in synergia_elements:
        element.set_double_attribute("aperture_radius", radius)

lattice_length = synergia_lattice.get_length()

reference_particle = synergia_lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
momentum = reference_particle.get_momentum()
brho = momentum / (synergia.foundation.pconstants.c / 1e9)
bunchlen_m = bunchlen_sec * beta * synergia.foundation.pconstants.c

emitx /= (beta * gamma) * 6.0
emity /= (beta * gamma) * 6.0

if myrank == 0:
    print "    lattice length                   :", lattice_length, "m"
    print "    brho                             :", brho, "T-m"
    print "    momentum                         :", momentum, "GeV/c"
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
freq = harmno * beta * synergia.foundation.pconstants.c / lattice_length
for element in synergia_elements:
    if element.get_type() == "rfcavity":
        element.set_double_attribute("volt", rf_voltage)
        element.set_double_attribute("freq", freq)
        if myrank == 0:
            print "    rfcavity                         :", element.get_name()
            print "    rf voltage                       :", 
            print element.get_double_attribute("volt"), "MV"
            print "    rf frequency                     :", 
            print element.get_double_attribute("freq") * 1e-6, "MHz"

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
t0 = time.time()
chef_beamline = chef_lattice.get_beamline()

#   Set up the tune control circuits
if myrank == 0:
    print
    print "....Set the tune control circuits"
horizontal_correctors = []
vertical_correctors = []
htest = []
index = 0
for element in synergia_elements:
    name = element.get_name()
    type = element.get_type()
    index += 1
    if name == "hqf1" or name == "hqf2":
        horizontal_correctors.append(element)
    if name == "hqd1" or name == "hqd2":
        vertical_correctors.append(element)

#   Set initial tunes.
if myrank == 0:
    print
    print "....Set initial tunes"
bare_tune = [opts.tuneh, opts.tunev]

#   Adjust Tunes (initial tunes)
if myrank == 0:
    print
    print "....Adjust Tunes"
tune_tolerance = 1.0e-8
adjuster_tune_step = 0.005

lattice_simulator.adjust_tunes(bare_tune[0], bare_tune[1], horizontal_correctors, vertical_correctors, tune_tolerance)

#   Store initial quad settings for initial tune
if myrank == 0:
    print
    print "....Store initial quad settings for initial tune"
initial_k1 = []
index = 0
for element in synergia_elements:
#lattice_simulator.get_lattice().get_elements():
    if element.get_type() == "quadrupole":
        initial_k1.append(element.get_double_attribute("k1"))
        index += 1

#   Output files for plotting
if opts.twiss and myrank == 0:
        twiss_file = ("twiss.txt")
        twiss_log = open(twiss_file, "w")
if opts.separatrix and myrank == 0:
        separatrix_file = ("separatrix.txt")
        separatrix_log = open(separatrix_file, "w")

#   Calculate target emittance from invariant emittance
#   (This is the area of a triangle tangent to a circle.)
if myrank == 0:
    print 
    print "....Calculate target emittance from invariant emittance"
safety_factor = 1.05
emittance_x = opts.norm_emit * numpy.pi / (beta * gamma)
emittance_x *= safety_factor 
emittance_x *= 3.0 * numpy.sqrt(3.0) / numpy.pi

#   Determine the target value of |g|
if myrank == 0:
    print
    print "....Determine the target value of |g|"
resonant_tune = opts.resonant_tune
delta = bare_tune[0] - resonant_tune
target_abs_g = numpy.abs(delta) / numpy.sqrt(2.0 * numpy.sqrt(3.0) * emittance_x)
target_phase_g = opts.phase_g * numpy.pi / 180.0
target_g = target_abs_g * (numpy.cos(target_phase_g) + 1j * numpy.sin(target_phase_g))
if myrank == 0:
    print "        resonant tune                :", resonant_tune
    print "        starting tune                :", bare_tune[0]
    print "        target |g|                   :", target_abs_g
    print "        target psi_g                 :", opts.phase_g, "degree"

#   Calculate current value of |g|
if myrank == 0:
    print
    print "....Calculate current value of g"
    print "........Calling resonance_sum"
    print

chef_beamline = lattice_simulator.get_chef_lattice().get_beamline()

for element in synergia_elements:
#lattice_simulator.get_lattice().get_elements():
    name = element.get_name()
    type = element.get_type()
    lattice_function = lattice_simulator.get_lattice_functions(element)

g = numpy.zeros([2], dtype=complex)
g = resonance_sum(lattice_simulator, brho)
total_g = g[0] + g[1]
abs_g = numpy.abs(total_g)

if myrank == 0:
    print
    print "        current value: g             :", total_g
    print "                      |g|            :", abs_g

#   Adjust harmonic sextupole strengths
if myrank == 0:
    print
    print "....Adjust harmonic sextupole strengths: Setting (integrated) sextupole strengths"
    print "........Change polarity, if desired"
if opts.flip:
    for element in chef_beamline:
        if element.Type() == "multipole":
            element.set_double_attribute("k2l", -1.0 * element.get_double_attribute("k2l"))

ratio = target_abs_g / abs_g
if myrank == 0:
    print 
    print "........Scale to target value of harmonic coupling: Setting (integrated) sextupole strengths"
    print "            ratio                    :", ratio
for element in synergia_elements:
    if element.get_type() == "multipole":
        old_strength = element.get_double_attribute("k2l")
        new_strength = old_strength * ratio
        element.set_double_attribute("k2l", new_strength)
        if myrank == 0:
            print 
            print "            type                     :", element.get_type()
            print "            name                     :", element.get_name()
            print "            strength                 :", 
            print old_strength * brho / 2.0,
            print "->", new_strength * brho / 2.0, "T/m"
            print "                                     : B''l =",
            print 2.0 * new_strength * brho / 2.0, "T/m"
g *= ratio

if myrank == 0:
    print 
    print "........Rotate, if desired"
denom = (numpy.conjugate(g[1]) * g[0]).imag
s = numpy.zeros([2], 'd')
s[0] = (numpy.conjugate(g[1]) * target_g).imag / denom;
s[1] = -(numpy.conjugate(g[0]) * target_g).imag / denom;
if myrank == 0:
    print "            g                        :", target_g
    print "            g1                       :", g[0]
    print "            g2                       :", g[1]
    print "            s1                       :", s[0]
    print "            s2                       :", s[1]
    print "            target_g                 :", numpy.abs(target_g)
    print "            g                        :", numpy.abs(g[0] + g[1])
    print "            target_phase_g           :", target_phase_g
    print "            delta                    :", delta
    print
    print "............Doing rotation"
for element in synergia_elements:
    if element.get_type() == "multipole":
        old_strength = element.get_double_attribute("k2l")

        if re.search('ddd_50', element.get_name()):
            new_strength = old_strength * s[0]
        elif re.search('ddd_20', element.get_name()):
            new_strength = old_strength * s[1]
        element.set_double_attribute("k2l", new_strength)
        if myrank == 0:
            print
            print "            type                     :", element.get_type()
            print "            name                     :", element.get_name()
            print "            strength                 :", 
            print old_strength * brho / 2.0,
            print "->", new_strength * brho / 2.0, "T/m"
            print "                                     : B''l =",
            print 2.0 * new_strength * brho / 2.0, "T/m"
g[0] *= s[0]
g[1] *= s[1]

#   Store final settings of harmonic sextupole circuits
final_k2l = []
index = 0
if myrank == 0:
    print
    print "........Store final settings of harmonic sextupole circuits"
for element in synergia_elements:
    if element.get_type() == "multipole":
        strength = element.get_double_attribute("k2l") / 2.0
        final_k2l.append(strength)        # 1/m^3
        index += 1
        if myrank == 0:
            print
            print "            type                     :", element.get_type()
            print "            name                     :", element.get_name()
            print "            strength (k2l)           :", final_k2l[index-1],
            print "1/m^3"
            print "                                     : B''l =",
            print strength * brho, "T/m"

sys.exit(1)



























###############################################################################
#   Set initial element strengths for beam matching
###############################################################################

###############################################################################
#  CHEF One Turn Map 
###############################################################################


###############################################################################
#   One Turn Map from Synergia
###############################################################################
map = linear_one_turn_map(lattice_simulator)
if myrank == 0:
    print "One turn map from synergia2.5 infrastructure"
    print numpy.array2string(map, max_line_width=200, precision=3)
    print "    det(M) :", numpy.linalg.det(map)

# checking eigen vector and values of map
#[l, v] = numpy.linalg.eig(map)
#print "l:", l
#print "v:", v
#if myrank == 0:
#    print "    eigenvalues:"
#    for z in l: 
#        print "    |z|:", abs(z), " z:", z, " tune:", numpy.log(z).imag/(2.0*numpy.pi)

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
    print "    xwidth                           :", numpy.sqrt(emitx * bx), "m"
    print "    ywidth                           :", numpy.sqrt(emity * by), "m"

if myrank == 0:
    print
    print "Begin generating bunch..."
#bunch = synergia.optics.generate_matched_bunch_transverse(lattice_simulator, 
#                emitx, emity, stdz, dpop, real_particles, macro_particles,
#                seed=seed)
bunch = synergia.optics.generate_matched_bunch(lattice_simulator,
                numpy.sqrt(emitx * beta_twiss[0]), numpy.sqrt(emity * beta_twiss[1]),
                bunchlen_m, real_particles, macro_particles, seed=seed)

# apply offset to bunch
particles = bunch.get_local_particles()

particles[:,0] = particles[:,0] + x_offset
particles[:,2] = particles[:,2] + y_offset
particles[:,4] = particles[:,4] + z_offset

if myrank == 0:
    print
    print "    expected stdx:", numpy.sqrt(emitx * bx), "(m) generated (on rank 0):", numpy.std(particles[:,0]), "(m)"
    print "    expected stdy:", numpy.sqrt(emity * by), "(m) generated (on rank 0):", numpy.std(particles[:,2]), "(m)"
    print "    expected stdz:", bunchlen_m, "   (m) generated (on rank 0):", numpy.std(particles[:,4]), "   (m)"

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
#propagator.propagate(bunch, num_turns, diagnostics_actions, ramp_actions, 
#                verbose)
t1 = time.time()

if myrank == 0:
    print
    print "propagate time =", t1 - t0
