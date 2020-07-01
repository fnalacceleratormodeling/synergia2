#!/usr/bin/env python

import sys, os
import numpy as np
import synergia

p_ke = 0.8 # 800 MeV injection energy
real_particles = 2.0e11
print "bunch charge (real_particles): ", real_particles
harmon = 50
lattice_length = 288.0
print "lattice length: ", lattice_length
p_energy = p_ke + synergia.foundation.pconstants.mp
print "p_energy: ", p_energy
p_momentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.mp, p_energy)
p_beta = p_momentum.get_beta()
p_gamma = p_momentum.get_gamma()
print "p_beta: ", p_beta
print "p_gamma: ", p_gamma

electron_ke = 0.00001  # 10 KeV electron beam
e_momentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.me)
e_momentum.set_kinetic_energy(electron_ke)
e_beta = e_momentum.get_beta()
print "e_beta: ", e_beta

bunch_length = lattice_length/harmon
print "bunch_length: ", bunch_length
bunch_current = synergia.foundation.pconstants.e * real_particles*p_beta*synergia.foundation.pconstants.c/bunch_length
print "bunch_current: ", bunch_current

longrms = 0.5
gaussian_current = synergia.foundation.pconstants.e * real_particles*p_beta*synergia.foundation.pconstants.c/(longrms*np.sqrt(2.0*np.pi))
print "gaussian_current: ", gaussian_current

bunch_current = gaussian_current

#sc_factor = 2.0 * bunch_current * lattice_length * synergia.foundation.pconstants.rp/(synergia.foundation.pconstants.e * synergia.foundation.pconstants.c * p_beta**3 * p_gamma**3)
#print "sc_factor: ", sc_factor

#elens_JL = sc_factor * p_beta**3 * p_gamma**3 * e_beta * p_beta**2 * p_gamma/(1.0 + e_beta*p_beta)

#elens_JL = bunch_current * lattice_length /(p_beta**3 * p_gamma**3) * e_beta * p_beta**2 * p_gamma/(1.0 + e_beta*p_beta)

elens_JL = bunch_current*lattice_length * e_beta /((1.0 + e_beta*p_beta)*p_beta*p_gamma**2)

print "peak elens_JL for full compensation: ", elens_JL

print "peak elens in 12 lens for full compensation: ", elens_JL/12.0

print "current for 6 lens each 2m long: ", elens_JL/(2.0*6)
print "current for 24 lens each 2m long: ", elens_JL/(2.0*24)
