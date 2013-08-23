#!/usr/bin/env python
import synergia
import numpy as np

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

lattice = synergia.lattice.Mad8_reader().get_lattice("model", "foborodobo32.lat")

# set frequency of RF cavities
lattice_length = lattice.get_length()          
reference_particle = lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()

# set rf cavity frequency based on length and harmonic number

# harmonic number of this lattice is 128
harmno = 128

# the stable frequency is harmno * beta * c/ring_length
freq = harmno * beta * synergia.foundation.pconstants.c/lattice_length
rfwavelen = lattice_length/harmno

print "rf frequency: ", freq

for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        # if the voltage was not set in the lattice file, you
        # would set it here.
        # elem.set_double_attribute("volt", rf_voltage)
        elem.set_double_attribute("freq", freq)
        # this lattice is above transition.  lag is specified as a fraction
        # of 2.0*pi
        elem.set_double_attribute("lag", 0.5)

synergia.lattice.xml_save_lattice(lattice, "foborodobo32_lattice.xml")

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)

map = synergia.optics.one_turn_map.linear_one_turn_map(lattice_simulator)

[ax, bx, qx] = map2twiss(map[0:2,0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6,4:6])

print "Lattice parameters (assuming uncoupled map)"
print "alpha_x: ", ax, " alpha_y: ", ay
print "beta_x: ", bx, " beta_y: ", by
print "q_x: ", qx, " q_y: ", qy
print "beta_z: ", bz

print "one turn map"
print np.array2string(map, max_line_width=300,precision=16)
