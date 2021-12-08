#!/usr/bin/env python3
import sys, os
import numpy as np
import synergia

def main():

    momentum = 10.0
    four_momentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.mp)
    four_momentum.set_momentum(momentum)
    refpart = synergia.foundation.Reference_particle(1, four_momentum)

    macroparticles = 100
    realparticles = 5e10
    num_bunch = 1
    bunch_spacing = 1.0

    commxx = synergia.utils.Commxx()

    lattice = synergia.lattice.Lattice("foo")
    d = synergia.lattice.Lattice_element("drift", "d")
    d.set_double_attribute("l", 10.0)
    lattice.append(d)
    print('Lattice has ', len(lattice.get_elements()), ' elements, length: ', lattice.get_length())

    print('the lattice:')
    print(lattice)
    print()
    
    sim = synergia.simulation.Bunch_simulator.create_bunch_train_simulator(
        refpart, macroparticles, realparticles, num_bunch, bunch_spacing)

    bunch = sim.get_bunch(0, 0)

    print('bunch has ', bunch.get_total_num(), ' particles')
    bunch.checkout_particles()
    particles = bunch.get_particles_numpy()
    localnum = bunch.get_local_num()
    print('particles.shape: ', particles.shape)
    print('num local particles: ', localnum)
    particles[:, 0:6] = 0.0

    bunch.checkin_particles()

    sim.reg_diag_per_step(synergia.bunch.Diagnostics_particles("foo.h5"), bunch_idx=0)

    sc_ops = synergia.collective.Dummy_CO_options()
    stepper = synergia.simulation.Split_operator_stepper(sc_ops, 1)

    print('lattice length: ', lattice.get_length())
    propagator = synergia.simulation.Propagator(lattice, stepper)

    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.DINFO)

    # check particles before propagate
    bunch.checkout_particles()
    particles = bunch.get_particles_numpy()
    print('first 2 particles')
    print(particles[:2,:])
    print()
    
    propagator.propagate(sim, simlog, 1)

    # check particles after propagate
    bunch.checkout_particles()
    particles = bunch.get_particles_numpy()
    print('first 2 particles')
    print(particles[:2,:])
    print()

if __name__ == "__main__":
    main()
