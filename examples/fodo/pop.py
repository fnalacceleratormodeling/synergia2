#!/usr/bin/env python3

from mpi4py import MPI
import numpy as np
import synergia


def print_statistics(bunch):

    parts = bunch.get_particles_numpy()
    print(parts.shape,  ", ", parts.size )
    print("shape: {0}, {1}".format(parts.shape[0], parts.shape[1]))

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch, mean)
    print("mean = {}".format(mean))
    print("std = {}".format(std))

def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            ref_part, 16384, 2.94e10)

    bunch = sim.get_bunch()

    bunch_means = np.zeros(6, dtype='d')
    # numerical values from cxx_covariance_matrix.xml in cxx_benchmark
    bunch_covariances = np.array(
        [[3.0509743977035345e-05, 2.2014134466660509e-06, 0, 0, 0, 0],
         [2.2014134466660509e-06, 1.9161816525115869e-07, 0, 0, 0, 0],
         [0, 0, 7.5506914064526925e-06, -6.6846812465678249e-07, 0, 0],
         [0, 0, -6.6846812465678249e-07, 1.9161816525115867e-07, 0, 0],
         [0, 0, 0, 0, 0.00016427607645871527, 0],
         [0, 0, 0, 0, 0, 1e-08]])
    
    dist = synergia.foundation.PCG_random_distribution(1234567, synergia.utils.Commxx())
    #dist = synergia.foundation.Random_distribution(1234567, synergia.utils.Commxx())
    synergia.bunch.populate_6d( dist, 
        bunch, 
        bunch_means,
        bunch_covariances)

    return sim

def run():

    screen = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.DEBUG)

    fourmomentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.mp)
    fourmomentum.set_momentum(3.0)
    refpart = synergia.foundation.Reference_particle(1, fourmomentum)
    print('before create_simulator')
    sim = create_simulator(refpart)
    print('after create_simulator')

    print('bunch statistics')
    print_statistics(sim.get_bunch())

def main():
    run()

main()
