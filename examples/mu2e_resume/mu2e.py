#!/usr/bin/env python
import sys
import os
import time
import synergia
import numpy
import re
import math
from synergia.optics.one_turn_map import linear_one_turn_map
from ramp_modules import Ramp_actions
from mu2e_options import opts
#from basic_toolkit import *
#from mxyzptlk import *
#from beamline import *
#from physics_toolkit import *

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

def get_resonance_sum(lattice_simulator, brho):
    imaginary_number = 1j
    factor = imaginary_number / ((6.0 * numpy.sqrt(2.0)) * (4.0 * numpy.pi) * brho)

    synergia_elements = lattice_simulator.get_lattice().get_elements()
    num_elements = len(synergia_elements)
    last_element = synergia_elements[num_elements-1]
    last_element_lattice_functions = lattice_simulator.get_lattice_functions(last_element)
    nu_x = last_element_lattice_functions.psi_x / (2.0 * numpy.pi)
    circumference = last_element_lattice_functions.arc_length
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
            #if myrank == 0:
            #    print
            #    print "        index                        :", index
            #    print "        type                         :", 
            #    print element.get_type()
            #    print "        name                         :", 
            #    print element.get_name()
            #    print "        strength                     :", 
            #    print element.get_double_attribute("k2l"), "1/m^2"
            #    print "        strength in chef unit        :",
            #    print element.get_double_attribute("k2l") * brho / 2.0, "T/m"
            #    print "        theta                        :", theta
            #    print "        beta_x                       :", beta_x
            #    print "        psi_x                        :", psi_x
            #    print "        phase                        :", phase
            increment = beta32 * brho \
                    * element.get_double_attribute("k2l") \
                    * (numpy.cos(phase) - numpy.sin(phase) * imaginary_number)
            #in debuncher.cc (ext_sim_e.cc)
            #std::complex<double> increment
            #=   std::complex<double>( beta32, 0.0 )
            #    * std::complex<double>( 2.0 * (*it)->Strength(), 0.0 )
            #    * std::complex<double>( cos(phase), - sin(phase) );
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
#   calculate lattice settings for third integer resonant extraction
###############################################################################
def calculate_lattice_settings():
    if myrank == 0:
        print
        print "Calculate and set lattice settings for the third integer resonant extraction"
    #   Set up the tune control circuits
    if myrank == 0:
        print
        print "....Set the tune control circuits"
    horizontal_correctors = []
    vertical_correctors = []
    index = 0 
    for element in synergia_elements:
        name = element.get_name()
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
    lattice_simulator.adjust_tunes(bare_tune[0], bare_tune[1], 
            horizontal_correctors, vertical_correctors, tune_tolerance)

    #   Store initial quad settings for initial tune
    if myrank == 0:
        print
        print "....Store initial quad settings for initial tune"
    initial_k1 = []
    index = 0
    for element in synergia_elements:
        if element.get_type() == "quadrupole":
            initial_k1.append(element.get_double_attribute("k1"))     # 1/m
            index += 1
            #if myrank == 0:
            #    print "       ", element.get_name(), initial_k1[index - 1], 
            #    print "1/m"

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
    target_abs_g = numpy.abs(delta) / numpy.sqrt(2.0 * numpy.sqrt(3.0) \
            * emittance_x)
    target_phase_g = opts.phase_g * numpy.pi / 180.0
    target_g = target_abs_g * (numpy.cos(target_phase_g) + 1j * \
            numpy.sin(target_phase_g))
    if myrank == 0:
        print "        resonant tune                :", resonant_tune
        print "        starting tune                :", bare_tune[0]
        print "        target |g|                   :", target_abs_g
        print "        target psi_g                 :", opts.phase_g, "degree"

    #   Calculate current value of |g|
    if myrank == 0:
        print
        print "....Calculate current value of g"
        print

    #~for element in synergia_elements:
    #~    name = element.get_name()
    #~    type = element.get_type()
    #~    lattice_function = lattice_simulator.get_lattice_functions(element)

    g = numpy.zeros([2], dtype=complex)
    g = get_resonance_sum(lattice_simulator, brho)
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
        for element in synergia_elements:
            if element.get_type() == "multipole":
                element.set_double_attribute("k2l", 
                        -1.0 * element.get_double_attribute("k2l"))

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
                print "            type                     :", 
                print element.get_type()
                print "            name                     :", 
                print element.get_name()
                print "            strength                 :",
                print new_strength, "1/m^2"
                print "            strength in chef unit    :",
                print old_strength * brho / 2.0,
                print "->", new_strength * brho / 2.0, "T/m"
                print "                                     : B''l =",
                print new_strength * brho, "T/m"
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
                print "            type                     :", 
                print element.get_type()
                print "            name                     :", 
                print element.get_name()
                print "            strength                 :",
                print new_strength, "1/m^2"
                print "            strength in chef unit    :",
                print old_strength * brho / 2.0,
                print "->", new_strength * brho / 2.0, "T/m"
                print "                                     : B''l =",
                print new_strength * brho, "T/m"
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
            strength = element.get_double_attribute("k2l")
            final_k2l.append(strength)        # 1/m^2
            index += 1
            if myrank == 0:
                print
                print "            type                     :", 
                print element.get_type()
                print "            name                     :", 
                print element.get_name()
                print "            strength (k2l)           :", 
                print final_k2l[index-1], "1/m^2"
                print "                                     : B''l =",
                print strength * brho, "T/m"

    #   Set resonant tune
    if myrank == 0:
        print
        print "....Set resonant tune"
    lattice_simulator.adjust_tunes(resonant_tune, bare_tune[1], 
            horizontal_correctors, vertical_correctors, tune_tolerance)

    #   Store final quad settings for resonant tune
    if myrank == 0:
        print
        print "....Store final quad settings for resonant tune"
    final_k1 = []
    index = 0
    for element in synergia_elements:
        if element.get_type() == "quadrupole":
            final_k1.append(element.get_double_attribute("k1"))         # 1/m
            index += 1
            #if myrank == 0:
            #    print "       ", element.get_name(), final_k1[index - 1], 
            #    print "1/m"

    #   Testing with final conditions
    if myrank == 0:
        print
        print "....Test: after setting final quad conditions"
    index = 0
    for element in synergia_elements:
        if element.get_type() == "quadrupole":
            element.set_double_attribute("k1", final_k1[index])
            index += 1
    lattice_simulator.update()
    map_final = linear_one_turn_map(lattice_simulator)
    qf = map2twiss(map_final[0:2, 0:2])[2]
    if myrank == 0:
        print "        qf                           :", qf
        print "        One Turn Map                 :"
        print numpy.array2string(map_final, max_line_width=200, precision=3)

    nu_h = lattice_simulator.get_horizontal_tune()
    nu_v = lattice_simulator.get_vertical_tune()
    fractional_tune = resonant_tune - int(resonant_tune)
    if myrank == 0:
        print
        print "        Fractional tunes: horizontal :", nu_h, 
        print "= resFracTune +", nu_h - fractional_tune
        print "                          vertical   :", nu_v

    #   Saving final lattice configurations
    if myrank == 0:
        print
        print "....Saving final lattice configurations"

    final_setting = ("/data/cspark/results/mu2e/%s/final_setting_tunes_%g_%g_phase_%g.xml" % (latt_dir, opts.tuneh, opts.tunev, opts.phase_g))
    synergia.lattice.xml_save_lattice(lattice_simulator.get_lattice(), 
            final_setting)

    #   Testing with initial conditions
    if myrank == 0:
        print
        print "....Test: after setting initial conditions"
    index = 0
    for element in synergia_elements:
        if element.get_type() == "quadrupole":
            element.set_double_attribute("k1", initial_k1[index])
            index += 1
    lattice_simulator.update()
    map_initial = linear_one_turn_map(lattice_simulator)
    qi = map2twiss(map_initial[0:2, 0:2])[2]
    if myrank == 0:
        print "        qi                           :", qi
        print "        One Turn Map                 :"
        print numpy.array2string(map_initial, max_line_width=200, precision=3)
    nu_h = lattice_simulator.get_horizontal_tune()
    nu_v = lattice_simulator.get_vertical_tune()
    fractional_tune = resonant_tune - int(resonant_tune)
    if myrank == 0:
        print
        print "        Fractional tunes: horizontal :", nu_h, 
        print "= resFracTune +", nu_h - fractional_tune
        print "                          vertical   :", nu_v

    #   Saving initial lattice configurations
    if myrank == 0:
        print
        print "....Saving initial lattice configuations"

    initial_setting = ("/data/cspark/results/mu2e/%s/initial_setting_tunes_%g_%g_phase_%g.xml" % (latt_dir, opts.tuneh, opts.tunev, opts.phase_g))
    synergia.lattice.xml_save_lattice(lattice_simulator.get_lattice(),
            initial_setting)

    return (initial_k1, final_k1, final_k2l)

###############################################################################
#   load pre-calculated lattice settings for third integer resonant extraction
###############################################################################
def load_lattice_settings():
    initial_setting = ("/data/cspark/results/mu2e/%s/initial_setting_tunes_%g_%g_phase_%g.xml" % (latt_dir, opts.tuneh, opts.tunev, opts.phase_g))
    initial_synergia_lattice = synergia.lattice.Lattice()
    synergia.lattice.xml_load_lattice(initial_synergia_lattice, initial_setting)
    initial_synergia_elements = initial_synergia_lattice.get_elements()
    #   Load initial quad settings
    if myrank == 0:
        print
        print "....Load initial quad settings"
    initial_k1 = []
    index = 0
    for element in initial_synergia_elements:
        if element.get_type() == "quadrupole":
            initial_k1.append(element.get_double_attribute("k1"))       # 1/m
            index += 1
            #if myrank == 0:
            #    print "       ", element.get_name(), initial_k1[index - 1], 
            #    print "1/m"

    final_setting = ("/data/cspark/results/mu2e/%s/final_setting_tunes_%g_%g_phase_%g.xml" % (latt_dir, opts.tuneh, opts.tunev, opts.phase_g))
    final_synergia_lattice = synergia.lattice.Lattice()
    synergia.lattice.xml_load_lattice(final_synergia_lattice, final_setting)
    final_synergia_elements = final_synergia_lattice.get_elements()
    #   Load final magnet settings
    if myrank == 0:
        print
        print "....Load final quad settings"
    final_k1 = []
    final_k2l = []
    index = 0
    index2 = 0
    for element in final_synergia_elements:
        if element.get_type() == "quadrupole":
            final_k1.append(element.get_double_attribute("k1"))         # 1/m
            index += 1
            #if myrank == 0:
            #    print "       ", element.get_name(), final_k1[index - 1], 
            #    print "1/m"
        if element.get_type() == "multipole":
            final_k2l.append(element.get_double_attribute("k2l"))       # 1/m^2
            index2 += 1
            if myrank == 0:
                print
                print "            type                     :",
                print element.get_type()
                print "            name                     :",
                print element.get_name()
                print "            strength (k2l)           :",
                print final_k2l[index2-1], "1/m^2"
                print "                                     : B''l =",
                print final_k2l[index2-1] * brho, "T/m"


    return (initial_k1, final_k1, final_k2l)

###############################################################################
#   calculate vertices for the star chamber
###############################################################################
def get_vertices(a, b):
    (u0, u1) = (2.2355, 0.0)
    (v0, v1) = (2.811, 2.811)
    (r0, r1) = (0.5755, 2.250 + 0.0435)
    (ang1, ang2) = (78.43, 66.86)
    d2r = numpy.pi / 180.0

    np = a + b + a - 1
    x = numpy.zeros(np,'d')
    y = numpy.zeros(np,'d')

    for i in range (0, a):
        j = i
        rot_ang = (i + 1) * ang1 / a
        x[j] = u0 + r0 * numpy.cos(rot_ang * d2r)
        y[j] = u1 + r0 * numpy.sin(rot_ang * d2r)

    for i in range(0, b):
        j = i + a
        rot_ang = (270 - (90 - ang1)) - (i + 1) * ang2 / b
        x[j] = v0 + r1 * numpy.cos(rot_ang * d2r)
        y[j] = v1 + r1 * numpy.sin(rot_ang * d2r)

    for i in range (0, a-1):
        j = i + a + b
        rot_ang = (90- ang1) + (i + 1) * ang1 / a
        x[j] = u1 + r0 * numpy.cos(rot_ang * d2r)
        y[j] = u0 + r0 * numpy.sin(rot_ang * d2r)

    return x * 0.0254, y * 0.0254

def get_all_vertices(x, y):
    a = 2.811 * 0.0254
    np_quad = x.size
    np = np_quad * 4 + 4
    xp = numpy.zeros(np,'d')
    yp = numpy.zeros(np,'d')
    (xp[0], yp[0])  = (a, 0.0)
    (xp[np_quad+1], yp[np_quad+1])  = (0.0, a)
    (xp[np_quad*2+2], yp[np_quad*2+2])  = (-a, 0.0)
    (xp[np_quad*3+3], yp[np_quad*3+3])  = (0.0, -a)
    for i in range(0, np_quad):
        index = i + np_quad * 0 + 1
        (xp[index], yp[index]) = (x[i], y[i])
    for i in range(0, np_quad):
        index = i + np_quad * 1 + 2
        (xp[index], yp[index]) = (-x[np_quad-i-1], y[np_quad-i-1])
    for i in range(0, np_quad):
        index = i + np_quad * 2 + 3
        (xp[index], yp[index]) = (-x[i], -y[i])
    for i in range(0, np_quad):
        index = i + np_quad * 3 + 4
        (xp[index], yp[index]) = (x[np_quad-i-1], -y[np_quad-i-1])
    return xp, yp

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

seed = opts.seed
grid = [gridx, gridy, gridz]
num_steps = opts.num_steps
num_turns = opts.num_turns
map_order = opts.map_order

radius = opts.radius
emitx = opts.norm_emit
emity = opts.norm_emit
bunchlen_sec = opts.bunchlen * 1e-9
rf_voltage = opts.rf_voltage
if rf_voltage == 0:
    latt_dir = "lattice_cache_norf_checkpoint"
else:
    latt_dir = "lattice_cache_rf_checkpoint"

x_offset = opts.x_offset
y_offset = opts.y_offset
z_offset = opts.z_offset

kick = opts.kick
wire_x = opts.wire_x
wire_width = opts.wire_width
gap = opts.gap

solver = opts.spacecharge

myrank = synergia.utils.Commxx().get_rank()

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
    print "    bunch length                     :", bunchlen_sec * 1e9, "nsec"
    print "    starting horizontal tune         :", opts.tuneh
    print "    starting vertical tune           :", opts.tunev
    print "    Radius aperture                  :", radius, "m"
    print "    RF Voltage                       :", rf_voltage, "MV"
    print "    X offset                         :", x_offset
    print "    Y offset                         :", y_offset
    print "    Z_offset                         :", z_offset
    print "    Septum kick strength (mrad)      :", kick * 1e3
    print "    Septum wire position (cm)        :", wire_x * 1e2
    print "    Septum wire width (um)           :", wire_width * 1e6
    print "    Septum wire gap (cm)             :", gap * 1e2
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

orig_lattice = synergia.lattice.Mad8_reader().get_lattice("debunch0",
                "Debunch_modified_rf.lat")
orig_elements = orig_lattice.get_elements()
synergia_lattice = synergia.lattice.Lattice()

#~septum_line = synergia.lattice.Mad8_reader().get_lattice("septum_line",
#~                "Debunch_modified_rf.lat")
#~septum_lattice_simulator = synergia.simulation.Lattice_simulator(
#~                septum_line, map_order)
#~septum_map = linear_one_turn_map(septum_lattice_simulator)
#~if myrank == 0:
#~    print "Transfer map for septum line"
#~    print numpy.array2string(septum_map, max_line_width=200, precision=3)
#~    print "    det(M) :", numpy.linalg.det(septum_map)
#~sys.exit(1)

#~lambertson_line = synergia.lattice.Mad8_reader().get_lattice(
#~                "lambertson_line", "Debunch_modified_rf.lat")
#~lambertson_lattice_simulator = synergia.simulation.Lattice_simulator(
#~                lambertson_line, map_order)
#~lambertson_map = linear_one_turn_map(lambertson_lattice_simulator)
#~if myrank == 0:
#~    print "Transfer map for lambertson line"
#~    print numpy.array2string(lambertson_map, max_line_width=200, precision=3)
#~    print "    det(M) :", numpy.linalg.det(lambertson_map)
#~sys.exit(1)

#   Element setup: apertures, septum, and lambertson
if myrank == 0:
    print
    print "Initial Element Setup"

np1 = 4
np2 = 10
total_np = (2*np1+np2-1)*4+4
(x, y) = get_vertices(np1, np2)
(u, v) = get_all_vertices(x, y)
pamr2 = numpy.min(u*u+v*v)

# Options for extractor type: "chef_propagate", "chef_map", or "chef_mixed"
for element in orig_elements:
    name = element.get_name()
    type = element.get_type()
    if type == "quadrupole":
        element.set_string_attribute("extractor_type", "chef_map")
        element.set_string_attribute("aperture_type","polygon")
        element.set_double_attribute("the_number_of_vertices", total_np)
        element.set_double_attribute("min_radius2", pamr2)
        for i in range(0, total_np):
            element.set_double_attribute("pax%d" % (i+1), u[i])
            element.set_double_attribute("pay%d" % (i+1), v[i])
        synergia_lattice.append(element)
    elif name == "e_septum":
        # generate new element for e_septum
        septum = synergia.lattice.Lattice_element("e_septum", name)
        # extractor type: e_septum must use Chef_propatator
        septum.set_string_attribute("extractor_type", "chef_propagate")
        # aperture type: e_septum uses wire_elliptical_aperture
        septum.set_string_attribute("aperture_type", "wire_elliptical")
        septum.set_double_attribute("wire_elliptical_aperture_horizontal_radius", radius)
        septum.set_double_attribute("wire_elliptical_aperture_vertical_radius", radius)
        septum.set_double_attribute("wire_elliptical_aperture_wire_x", -wire_x)
        septum.set_double_attribute("wire_elliptical_aperture_wire_width", wire_width)
        septum.set_double_attribute("wire_elliptical_aperture_gap", gap)
        # settings for e_septum wire and kick
        septum.set_double_attribute("positive_strength", 0.0);
        septum.set_double_attribute("negative_strength", kick);
        septum.set_double_attribute("wire_position", wire_x);
        septum.set_double_attribute("wire_width", wire_width);
        septum.set_double_attribute("gap_size", gap);
        septum.set_string_attribute("force_diagnostics", "true")
        synergia_lattice.append(septum)
    elif name == "lambertson":
        # generate new element for lambertson
        lambertson = synergia.lattice.Lattice_element("lambertson", name)
        # extractor type: lambertson must use Chef_propatator
        lambertson.set_string_attribute("extractor_type", "chef_propagate")
        # aperture type: lambertson uses lambertson_aprture
        lambertson.set_string_attribute("aperture_type", "lambertson")
        lambertson.set_double_attribute("lambertson_aperture_radius", -wire_x)
        lambertson.set_string_attribute("force_diagnostics", "true")
        synergia_lattice.append(lambertson)
    elif type == "multipole":
        element.set_string_attribute("extractor_type", "chef_propagate")
        element.set_string_attribute("aperture_type","circular")
        element.set_double_attribute("circular_aperture_radius", radius)
        synergia_lattice.append(element)
    elif type == "sextupole":
        element.set_string_attribute("extractor_type", "chef_propagate")
        element.set_string_attribute("aperture_type","circular")
        element.set_double_attribute("circular_aperture_radius", radius)
        synergia_lattice.append(element)
    elif type == "sbend":
        element.set_string_attribute("extractor_type", "chef_propagate")
        element.set_string_attribute("aperture_type","circular")
        element.set_double_attribute("circular_aperture_radius", radius)
        synergia_lattice.append(element)
    else:
        element.set_string_attribute("extractor_type", "chef_propagate")
        element.set_string_attribute("aperture_type","circular")
        element.set_double_attribute("circular_aperture_radius", radius)
        synergia_lattice.append(element)

synergia_elements = synergia_lattice.get_elements()
synergia_lattice.set_reference_particle(orig_lattice.get_reference_particle())

#index = 0
#length = 0
#for element in synergia_elements:
#    index += 1
#    name = element.get_name()
#    type = element.get_type()
#    length += element.get_length()
#    if myrank == 0:
#        print index, name, type, length
#sys.exit(1)

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

resonant_tune = opts.resonant_tune
delta = opts.tuneh - resonant_tune

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
#   rf cavity
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
#   Set lattice_simulator
###############################################################################
if myrank == 0:
    print
    print "Set lattice_simulator"

lattice_simulator = synergia.simulation.Lattice_simulator(synergia_lattice, 
                map_order)

#if myrank == 0:
#    index = 0
#    for element in lattice_simulator.get_chef_lattice().get_beamline():
#        index += 1
#        if element.Type() != "marker":
#            print index, element.Length(), element.Name()
#sys.exit(1)

###############################################################################
#   Ramping Actions
###############################################################################
#   load/calculate initial and final magnet settings
t0 = time.time()
if opts.lattice_load == 1:
    [initial_k1, final_k1, final_k2l] = load_lattice_settings()
else:
    [initial_k1, final_k1, final_k2l] = calculate_lattice_settings()
t1 = time.time()
if myrank == 0:
    print
    print "....Lattice setting time =", t1 - t0

max_turns = opts.max_turns
ramp_turns = opts.rampturns
extraction_fraction = opts.extraction_fraction
#~Below is the old method to calculate the number of turns to extract the 
#~whole particles.
#~turns_to_extract = int(extraction_fraction * beta \
#~        * synergia.foundation.pconstants.c / lattice_length / 60.0 ) \
#~        + ramp_turns
turns_to_extract = int(opts.spill_time * 1e-3 * beta \
        * synergia.foundation.pconstants.c / lattice_length) + ramp_turns

if myrank == 0:
    print
    print "....Determine ramping turns"
    print "        sextupole ramping turns      :", ramp_turns
    print "        extraction fraction          :", extraction_fraction
    print "        turns to extract:            :", turns_to_extract
    print "        final turn                   :", num_turns

ramp_actions = Ramp_actions(ramp_turns, turns_to_extract, initial_k1,
                final_k1, final_k2l)

###############################################################################
#   Set initial element strengths for beam matching
###############################################################################
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
lattice_simulator.update()

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


###############################################################################
#   Generating Bunch
###############################################################################
[ax, bx, qx] = map2twiss(map[0:2, 0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6, 4:6])
alpha_twiss, beta_twiss = synergia.optics.get_alpha_beta(map)
if ax != alpha_twiss[0] or ay != alpha_twiss[1] or bx != beta_twiss[0] or by != beta_twiss[1]:
    print "Error:"    #TODO
    sys.out(1)
dpop = bunchlen_m / bz

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
#        emitx, emity, stdz, dpop, real_particles, macro_particles,
#        seed=seed)
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

#Bunch_simulator
bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

###############################################################################
#   Collective operator
###############################################################################
if myrank == 0:
    print
    print "Set collective operator"

if opts.comm_avoid:
    collective_comm = synergia.utils.Commxx(True)
else:
    collective_comm = synergia.utils.Commxx()

if solver == "3d" or solver == "3D":
    if myrank == 0:
        print "    using 3D Open Hockney space charge solver"
    space_charge = synergia.collective.Space_charge_3d_open_hockney(
                    collective_comm, grid, False)
    stepper = synergia.simulation.Split_operator_stepper(lattice_simulator,
                    space_charge, num_steps)
elif solver == "2d" or solver == "2D":
    if myrank == 0:
        print "    using 2D Open Hockney space charge solver"
    space_charge = synergia.collective.Space_charge_2d_open_hockney(
                    collective_comm, grid, opts.convert_state, False)
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
#~old method for diagnostics
#~diagnostics_actions = synergia.simulation.Standard_diagnostics_actions()

#diagnostics per step
for part in range(0, opts.step_tracks):
    bunch_simulator.add_per_step(synergia.bunch.Diagnostics_track(
                            "mu2e_step_track_%02d.h5" % part, part))
if opts.step_full2:
    bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2(
                            "mu2e_step_full2.h5"))
if opts.step_particles:
    bunch_simulator.add_per_step(synergia.bunch.Diagnostics_particles(
                            "mu2e_step_particles.h5"))
if opts.forced_diagnostics:
    bunch_simulator.add_per_forced_diagnostics_step(
                            synergia.bunch.Diagnostics_full2("mu2e_forced_full2.h5"))

# diagnostics per turn
for part in range(0, opts.turn_tracks):
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_track(
                            "mu2e_turn_track_%02d.h5" % part, part))
if opts.turn_full2:
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2(
                            "mu2e_turn_full2.h5"))

if opts.turn_particles:
    bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles(
                            "mu2e_turn_particles.h5"), opts.particles_period)

###############################################################################
#   Propagate and track particles
###############################################################################
if myrank == 0:
    print
    print "Begin propagate..."

if max_turns <= ramp_turns:
    ramp_max_turns = opts.max_turns
else:
    ramp_max_turns = ramp_turns


t0 = time.time()
propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(opts.checkpointperiod)
#propagator.set_checkpoint_with_xml(True)
propagator.propagate(bunch_simulator, ramp_actions, num_turns, ramp_max_turns,
                opts.verbosity)
t1 = time.time()

lattice_diagnostics = synergia.lattice.Lattice_diagnostics(synergia_lattice,
                "lattice_deposited_charge_ramp.h5", "deposited_charge")
lattice_diagnostics.update_and_write()

if max_turns > ramp_turns:
    del lattice_diagnostics
    if myrank == 0:
        os.system("mv log log_ramp")

    if myrank == 0:
        print
        print "Propagate sextupole ramping finished"
        print "Propagate time =", t1 - t0

    resume = synergia.simulation.Resume()
    content = resume.get_content()
    resume.propagate(True, max_turns-ramp_turns, False, opts.verbosity)

    t1 = time.time()

    lattice_diagnostics = synergia.lattice.Lattice_diagnostics(content.lattice,
                    "lattice_deposited_charge.h5", "deposited_charge")
    lattice_diagnostics.update_and_write()

if myrank == 0:
    print
    print "Finish propagate turns"
    print "Propagate time =", t1 - t0
