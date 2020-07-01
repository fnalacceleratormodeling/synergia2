#!/usr/bin/env python

import sys, os
import numpy as np
import synergia

lattice1 = synergia.lattice.Lattice(synergia.utils.read_lsexpr_file("adjusted_lattice.lsx"))
#for elem in lattice1.get_elements():
#    elem.set_string_attribute("extractor_type", "chef_map")

lattice2x = synergia.lattice.MadX_reader().get_lattice("machine", "adjusted_lattice.madx")


lattice2 = synergia.lattice.Lattice("lattice2")
lattice2.set_reference_particle(lattice2x.get_reference_particle())
for elem in lattice2x.get_elements():
#    elem.set_string_attribute("extractor_type", "chef_map")
    lattice2.append(elem)

macro_particles = 1024
real_particles = 1.0e10

commxx = synergia.utils.Commxx()

refpart1 = lattice1.get_reference_particle()
energy1 = refpart1.get_total_energy()
refpart2 = lattice2.get_reference_particle()
energy2 = refpart2.get_total_energy()

if energy1 != energy2:
    raise RuntimeError, "energies of two lattices do not agree! %g != %g, diff: %g"%(energy1, energy2, (energy1-energy2)/energy2)

energy = refpart1.get_total_energy()
momentum = refpart1.get_momentum()
gamma = refpart1.get_gamma()
beta = refpart1.get_beta()
print "beam energy: ", energy
print "beam momentum: ", momentum
print "beam gamma: ", gamma
print "beam beta: ", beta

stepper1 = synergia.simulation.Independent_stepper_elements(lattice1, 1, 1)
lattice_simulator1 = stepper1.get_lattice_simulator()
stepper2 = synergia.simulation.Independent_stepper_elements(lattice2, 1, 1)
lattice_simulator2 = stepper2.get_lattice_simulator()

lf1 = []
for elem in lattice1.get_elements():
    lf1.append(lattice_simulator1.get_lattice_functions(elem))

emitx = 1.00055865e-6
beta_x1 = lf1[-1].beta_x
beta_y1 = lf1[-1].beta_y

lf2 = []
for elem in lattice2.get_elements():
    lf2.append(lattice_simulator2.get_lattice_functions(elem))

emitx = 1.00055865e-6
beta_x1 = lf1[-1].beta_x
beta_y1 = lf1[-1].beta_y
beta_x2 = lf2[-1].beta_x
beta_y2 = lf2[-1].beta_y

print "initial beta_x1: ", beta_x1
print "initial beta_x2: ", beta_x2
print "initial beta_y1: ", beta_y1
print "initial beta_y2: ", beta_y2

if abs(beta_x1/beta_x2-1) > 4.0e-14:
    raise RuntimeError, "beta_x functions differ: %.17g != %.17g: delta = %g"%(beta_x1, beta_x2, (beta_x1/beta_x2 - 1.0))

if abs(beta_y1/beta_y2-1) > 4.0e-14:
    raise RuntimeError, "beta_y functions differ: %.17g != %.17g: delta = %g"%(beta_y1, beta_y2, (beta_y1/beta_y2 - 1.0))

stdx = np.sqrt(beta_x1 * emitx)
stdy = np.sqrt(beta_y1 * emitx)
zrms = 0.5

bunch1 = synergia.optics.generate_matched_bunch(lattice_simulator1, stdx, stdy, zrms, real_particles, macro_particles, rms_index=[0,2,4], periodic=False, seed=1234567)

diag_particles = synergia.bunch.Diagnostics_particles("init_particles.h5")
diag_particles.set_bunch(bunch1)
diag_particles.update_and_write()

bunch2 = synergia.bunch.Bunch(refpart2, macro_particles, real_particles, commxx)
lp2 = bunch2.get_local_particles()
lp1 = bunch1.get_local_particles()

lp2[:,:] = lp1[:,:]

bunch_simulator1 = synergia.simulation.Bunch_simulator(bunch1)
bunch_simulator2 = synergia.simulation.Bunch_simulator(bunch2)
bunch_simulator1.add_per_step(synergia.bunch.Diagnostics_particles("p1.h5"))
bunch_simulator2.add_per_step(synergia.bunch.Diagnostics_particles("p2.h5"))
bunch_simulator1.add_per_step(synergia.bunch.Diagnostics_bulk_track("tracks1.h5", macro_particles))
bunch_simulator2.add_per_step(synergia.bunch.Diagnostics_bulk_track("tracks2.h5", macro_particles))

propagator1 = synergia.simulation.Propagator(stepper1)
propagator2 = synergia.simulation.Propagator(stepper2)

propagator1.propagate(bunch_simulator1, 1, 1, 1)
propagator2.propagate(bunch_simulator2, 1, 1, 1)
