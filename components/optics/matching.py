#!/usr/bin/env python

import numpy

from one_turn_map import linear_one_turn_map
from mpi4py import MPI
from pybunch import Bunch, populate_6d
from pyfoundation import Random_distribution
from math import acos, sin, sqrt

def _get_correlation_matrix(map, stdx, stdy, stdz):
    evals, evect_matrix = numpy.linalg.eig(map)
    evects = []
    for i in range(0, 6):
        evects.append(evect_matrix[:, i])
    F = range(0, 3)
    remaining = range(5, -1, -1)
    for i in range(0, 3):
        # find complex conjugate among remaining eigenvectors
        first = remaining.pop()
        best = 1.0e30
        conj = -1
        for item in remaining:
            sum = evects[first] + evects[item]
            if abs(numpy.max(sum.imag)) < best:
                best = abs(numpy.max(sum.imag))
                conj = item
        if conj == -1:
            raise RuntimeError, "failed to find a conjugate pair in ha_match"
        remaining.remove(conj)
        tmp = numpy.outer(evects[first],
            numpy.conjugate(evects[first]))
        tmp += numpy.outer(evects[conj],
            numpy.conjugate(evects[conj]))
        F[i] = tmp.real

    S = numpy.zeros((3, 3), 'd')
    for i in range(0, 3):
        for j in range(0, 3):
            S[i, j] = F[j][i, i]

    Sinv = numpy.linalg.inv(S)

    C = numpy.zeros([6, 6], 'd')
    units = [1, 1, 1, 1, 1, 1] # jfa have to think about units!
    cd1 = stdx * units[0] * stdx * units[0]
    cd2 = stdy * units[1] * stdy * units[1]
    cd3 = stdz * units[2] * stdz * units[2]

    for i in range(0, 3):
        C += F[i] * (Sinv[i, 0] * cd1 + Sinv[i, 1] * cd2 + Sinv[i, 2] * cd3)
    return C

def get_alpha_beta(map):
    u = numpy.ones([6])
    mxx = map[0, 0]
    mxpxp = map[1, 1]
    mxxp = map[0, 1]
    cos_mu = (mxx + mxpxp) / 2.0
    mu = acos(cos_mu)
    # beta function is positive
    # use this to pick branch
    if mxxp / sin(mu) < 0:
        mu = 2 * pi - mu
    beta_x = mxxp / sin(mu) * u[1] / u[0]
    alpha_x = (mxx - mxpxp) / (2.0 * sin(mu))

    myy = map[2, 2]
    mypyp = map[3, 3]
    myyp = map[2, 3]
    cos_mu = (myy + mypyp) / 2.0
    mu = acos(cos_mu)
    # beta function is positive
    # use this to pick branch
    if myyp / sin(mu) < 0:
        mu = 2 * pi - mu
    beta_y = myyp / sin(mu) * u[3] / u[2]
    alpha_y = (myy - mypyp) / (2.0 * sin(mu))

    return [alpha_x, alpha_y], [beta_x, beta_y]

def match_transverse_twiss_emittance(emittance, alpha, beta):
    """Calculate input parameters for a matched beam of given width
    using Courant-Snyder (Twiss) parameters. Returns
        sigma4,r4
    where
        sigma4 = [sigmax, sigmaxp, sigmay, sigmayp]
    and
        r4 = [rxxp, ryyp]
    """
    sigma4 = [0, 0, 0, 0]
    r4 = [0, 0]
    for i in range(0, 2):
        gamma = (1 + alpha[i] ** 2) / beta[i]
        width = sqrt(beta[i] * emittance[i])
        width_prime = sqrt(gamma * emittance[i])
        sigma4[2 * i] = width
        sigma4[2 * i + 1] = width_prime
        r4[i] = -alpha[i] / sqrt(1 + alpha[i] ** 2)
    return sigma4, r4

def get_covariances(sigma, r):
    c = numpy.zeros((6, 6), 'd')
    for i in range(0, 6):
        c[i, i] = sigma[i] ** 2
    for (i, ri) in zip([0, 2, 4], r):
        c[i, i + 1] = ri * sigma[i] * sigma[i + 1]
    return c

def generate_matched_bunch(lattice_simulator, stdx, stdy, stdz,
                           num_real_particles, num_macro_particles, seed=0,
                           comm=None):

    map = linear_one_turn_map(lattice_simulator)
    correlation_matrix = _get_correlation_matrix(map, stdx, stdy, stdz)
    if comm == None:
        comm = MPI.COMM_WORLD
    bunch = Bunch(lattice_simulator.get_lattice().get_reference_particle(),
                  num_macro_particles, num_real_particles, comm);

def generate_matched_bunch_transverse(lattice_simulator, emit_x, emit_y, rms_z, dpop,
                           num_real_particles, num_macro_particles, seed=0,
                           comm=None):

    map = linear_one_turn_map(lattice_simulator)
    alpha, beta = get_alpha_beta(map)
    sigma4, r4 = match_transverse_twiss_emittance([emit_x, emit_y], alpha, beta)
    sigma = range(0, 6)
    sigma[0:4] = sigma4
    beta = lattice_simulator.get_lattice().get_reference_particle().get_beta()
    p = lattice_simulator.get_lattice().get_reference_particle().get_momentum()
    sigma[4] = rms_z / beta
    sigma[5] = dpop * p
    r = [r4[0], r4[1], 1.0]
    covariance_matrix = get_covariances(sigma, r)
    means = numpy.zeros((6,), 'd')
    if comm == None:
        comm = MPI.COMM_WORLD
    bunch = Bunch(lattice_simulator.get_lattice().get_reference_particle(),
                  num_macro_particles, num_real_particles, comm)
    dist = Random_distribution(seed, comm)
    populate_6d(dist, bunch, means, covariance_matrix)
    return bunch
