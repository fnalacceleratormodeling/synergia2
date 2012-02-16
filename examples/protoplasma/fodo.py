#!/usr/bin/env python
import synergia
import sys
import numpy
import mpi4py.MPI as MPI
from fodo_options import opts

rank = MPI.COMM_WORLD.rank

#~ For twiss parameter calculations
#lattice_sextant=synergia.lattice.Mad8_reader().get_lattice("lo_beta_sextant", "tevatron.lat")
#sextant_elements = lattice_sextant.get_elements()
#sextant_simulator = synergia.simulation.Lattice_simulator(lattice_sextant, opts.map_order)
#if myrank == 0:
#    twiss_file = ("twiss.txt")
#    twiss_log = open(twiss_file, "w")
#index = 0
#length = 0.0
#for element in sextant_elements:
#    lattice_functions = sextant_simulator.get_lattice_functions(element)
#    type = element.get_type()
#    name = element.get_name()
#    length += element.get_length()
#    beta_x = lattice_functions.beta_x
#    beta_y = lattice_functions.beta_y
#    alpha_x = lattice_functions.alpha_x
#    D_x = lattice_functions.D_x
#    if myrank == 0:
#        twiss_log.write("%3d %12s %10s %10.5f %10.5f %10.5f %10.5f %10.5f\n" % (
#                    index, type, name, length, beta_x, beta_y, alpha_x, D_x))
#        twiss_log.flush()
#    index += 1
#if myrank == 0:
#    twiss_log.close()

lattice_Tev = synergia.lattice.Mad8_reader().get_lattice("injection", "tevatron.lat")

lattice_simulator_Tev = synergia.simulation.Lattice_simulator(lattice_Tev, opts.map_order)

bunch = synergia.optics.generate_matched_bunch_transverse(
              lattice_simulator_Tev, opts.emit, opts.emit, opts.stdz, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)

lattice_length = lattice_Tev.get_length()
reference_particle = lattice_Tev.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
p_ref = reference_particle.get_momentum() 
mass = reference_particle.get_four_momentum().get_mass()
brho = p_ref / (synergia.foundation.pconstants.c / 1e9)
if rank == 0:
    print
    print "    Injection line(hi_beta_sextant) information"
    print "    lattice length                   :", lattice_length, "m"
    print "    brho                             :", brho, "T-m"
    print "    momentum                         :", p_ref, "GeV/c"
    print "    energy                           :", energy, "GeV"
    print "    beta                             :", beta
    print "    gamma                            :", gamma

particles = bunch.get_local_particles()
local_num = bunch.get_local_num()
c = synergia.foundation.pconstants.c
half_num = int(local_num/2)
for part in range(0,local_num):
    p = (1.0 + particles[part,5]) * p_ref
    eoc = numpy.sqrt(p * p + mass * mass)
    if part >= half_num:
        epoc = eoc + opts.deltae
    else:
        epoc = eoc - opts.deltae
    p_prime = numpy.sqrt(epoc * epoc - mass * mass)
    particles[part,5] = (p_prime - p_ref) / p_ref

lattice_line = synergia.lattice.Mad8_reader().get_lattice("plasma", "tevatron.lat")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice_line, opts.map_order)

lattice_length = lattice_line.get_length()
reference_particle = lattice_line.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
p_ref = reference_particle.get_momentum()
brho = p_ref / (synergia.foundation.pconstants.c / 1e9)

if rank == 0:
    print
    print "    Plasma line(lo_beta_straight) information"
    print "    lattice length                   :", lattice_length, "m"
    print "    brho                             :", brho, "T-m"
    print "    momentum                         :", p_ref, "GeV/c"
    print "    energy                           :", energy, "GeV"
    print "    beta                             :", beta
    print "    gamma                            :", gamma

stepper = synergia.simulation.Independent_stepper_elements(
                            lattice_simulator, opts.steps)
diagnostics_step = synergia.bunch.Diagnostics_full2(bunch, "full2.h5")
diagnostics_turn = synergia.bunch.Diagnostics_particles(bunch, "particles.h5")
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, opts.turns,
                     diagnostics_step, diagnostics_turn,
                     opts.verbose)
