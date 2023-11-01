#!/usr/bin/env python
import sys
import os
import numpy as np
import synergia

# Example demonstrating how to emit diagnostics on a bunch at
# an arbitrary time

# Lorentz variables
# betagamma**2 + 1 = gamma**2
# (3/4)**2 + (4/4)**2 = (5/4)**2
betagamma = 3/4
gamma = 5/4
beta = betagamma/gamma

macroparticles = 10000
realparticles = 5e10

#######################################################

def print_statistics(bunch, fout=sys.stdout):

    parts = bunch.get_particles_numpy()
    print(parts.shape,  ", ", parts.size , file=fout)
    print("shape: {0}, {1}".format(parts.shape[0], parts.shape[1]))

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch, mean)
    print("mean = {}".format(mean), file=fout)
    print("std = {}".format(std), file=fout)

#######################################################

# get a covariance matrix
def get_covar_matrix():
    beta_x = 5.0
    alpha_x = 0.0
    beta_y = 30.0
    alpha_y = 0.0
    beta_z = 1000.0

    emitxy = 1.0e-6
    var_x = emitxy*beta_x
    var_xp = emitxy/beta_x
    var_y = emitxy*beta_y
    var_yp = emitxy/beta_y
    covxxp = 0.0 # apha = 0
    covyyp = 0.0 # alpha = 0
    var_z = 0.25
    var_cdt = var_z/beta**2
    var_dpop = var_cdt/(beta_z**2)

    M = np.array([ [var_x, 0, 0, 0, 0, 0],
                  [0, var_xp, 0, 0, 0, 0],
                  [0, 0, var_y, 0, 0, 0],
                  [0, 0, 0, var_yp, 0, 0],
                  [0, 0, 0, 0, var_cdt, 0],
                  [0, 0, 0, 0, 0, var_dpop]])
    return M

#######################################################


# Generate the test bunch
def get_bunch():
    # should use a real particle mass because openPMD might get confused
    commxx = synergia.utils.Commxx()
    refpart =  synergia.foundation.Reference_particle(1,
                                                      synergia.foundation.pconstants.mp,
                                                      gamma*synergia.foundation.pconstants.mp)
    bunch = synergia.bunch.Bunch(refpart, macroparticles, realparticles, commxx)
    dist = synergia.foundation.Random_distribution(1234579, commxx.get_rank()) 
    covars = get_covar_matrix()
    means = np.zeros(6, dtype='d')
    synergia.bunch.populate_6d(dist, bunch, means, covars)
    return bunch

#######################################################

def main():
    bunch = get_bunch()

    diag = synergia.bunch.Diagnostics_particles('particles.h5')
    handler,id =  bunch.add_diagnostics(diag)
    handler.update_and_write()
    print_statistics(bunch)

#######################################################

if __name__ == "__main__":
    main()

