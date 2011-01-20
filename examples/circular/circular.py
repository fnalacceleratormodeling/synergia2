#!/usr/bin/env python
import synergia
from circular_options import opts
import numpy as np
from synergia.optics.one_turn_map import linear_one_turn_map

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

lattice = synergia.lattice.Mad8_reader().get_lattice("model", "foborodobo_s.lat")
lattice_length = lattice.get_length()
print "lattice length: ", lattice_length
reference_particle = lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()

print "reference particle energy: ", energy
print "reference particle beta: ", beta
print "reference particle gamma: ", gamma

# set rf cavity frequency
# harmno * beta * c/ring_length
freq = opts.harmno * beta * synergia.foundation.pconstants.c/lattice_length
print "RF freq: ", freq

# Don't need that?  Lattice has voltage set
# rf cavity voltage, 
for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", opts.rf_voltage)
        elem.set_double_attribute("freq", freq)
        #elem.set_double_attribute("lag", 0.5)

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)

map = linear_one_turn_map(lattice_simulator)


print "Linear one turn map"
print np.array2string(map,max_line_width=200)

l,v = np.linalg.eig(map)
print "eigenvalues of one turn map: ", l
print "absolute values of eigenvalues (should all be 1): ", abs(l)
print "fractional tunes from eigenvalues: ", np.log(l).imag/(2.0*np.pi)


[[ax,ay],[bx,by]] = synergia.optics.get_alpha_beta(map)
print "Lattice functions assuming uncoupled map:"
print "alpha x: ", ax
print "alpha y: ", ay
print "beta x: ", bx
print "beta y: ", by

[az, bz, qz] = map2twiss(map[4:6,4:6])
print "alpha z (better be small): ", az
print "beta z: ", bz


no_op = synergia.simulation.Dummy_collective_operator("stub")
stepper = synergia.simulation.Split_operator_stepper(
                            lattice_simulator, no_op, opts.num_steps)

emit = opts.norm_emit
print "generating particles with transverse emittance: ", emit

covar = synergia.optics.matching._get_correlation_matrix(map, np.sqrt(emit*bx), np.sqrt(emit*by), opts.stdz)
print "covariance matrix"
print np.array2string(covar,max_line_width=200)


bunch = synergia.optics.generate_matched_bunch(lattice_simulator,
                                               np.sqrt(emit*bx), np.sqrt(emit*by), opts.stdz,
                                               opts.num_real_particles,
                                               opts.num_macro_particles,
                                               seed=opts.seed)
                                               
# get particle bunch for examination
particles = bunch.get_local_particles()

# apply offset to bunch
particles[:,0] = particles[:,0]+opts.x_offset
particles[:,2] = particles[:,2]+opts.y_offset
particles[:,4] = particles[:,4]+opts.z_offset

print "expect std_x: ", np.sqrt(emit*bx)
print "generated std_x: ", np.std(particles[:,0])
print "expect std_y: ", np.sqrt(emit*by)
print "generated std_y: ", np.std(particles[:,2]);
print "expected std_z: ", opts.stdz
print "generated std_z: ", np.std(particles[:,4])
print "expected std(dpop): ", opts.stdz/bz
print "generated std(dpop): ", np.std(particles[:,5])

diagnostics_writer_step = synergia.bunch.Diagnostics_writer("circular_full2.h5",
                                                            synergia.bunch.Diagnostics_full2())
diagnostics_writer_turn = synergia.bunch.Diagnostics_writer("circular_particles.h5",
                                                            synergia.bunch.Diagnostics_particles())
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch, opts.num_turns, diagnostics_writer_step, diagnostics_writer_turn, opts.verbose)
