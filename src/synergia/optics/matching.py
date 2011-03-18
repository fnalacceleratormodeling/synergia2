#!/usr/bin/env python
import sys
import numpy

from one_turn_map import linear_one_turn_map
from mpi4py import MPI
from synergia.bunch import Bunch, populate_6d, populate_transverse_gaussian
from synergia.foundation import Random_distribution
from math import acos, pi, sin, sqrt

def _get_correlation_matrix(map, rms_x, rms_y, rms_z, beta):
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
        # F[i] is effectively 2*e[i] cross e^H[i].

    # The correlation matrix is a linear combination of F[i] with
    # appropriate coefficients such that the diagonal elements C[i,i] i=(0,2,4)
    # come out to be the desired 2nd moments.
    S = numpy.zeros((3, 3), 'd')
    for i in range(0, 3):
        for j in range(0, 3):
            S[i, j] = F[j][2 * i, 2 * i]

    Sinv = numpy.linalg.inv(S)

    C = numpy.zeros([6, 6], 'd')
    cd1 = rms_x ** 2
    cd2 = rms_y ** 2
    rms_cdt = rms_z / beta
    cd3 = rms_cdt ** 2

    for i in range(0, 3):
        C += F[i] * (Sinv[i, 0] * cd1 + Sinv[i, 1] * cd2 + Sinv[i, 2] * cd3)
    return C

def get_alpha_beta(map):
    u = numpy.ones([6])
    mxx = map[0, 0]
    mxpxp = map[1, 1]
    mxxp = map[0, 1]
    cos_mu = (mxx + mxpxp) / 2.0
    try:
        mu = acos(cos_mu)
    except ValueError, e:
        sys.stderr.write('failed to take the acos of %g\n' % cos_mu)
        raise RuntimeError("get_alpha_beta: unstable map:\n" +
                           numpy.array2string(map) + '\n')
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

def generate_matched_bunch(lattice_simulator, rms_x, rms_y, rms_z,
                           num_real_particles, num_macro_particles, seed=0,
                           comm=None):

    map = linear_one_turn_map(lattice_simulator)
    beta = lattice_simulator.get_lattice().get_reference_particle().get_beta()
    correlation_matrix = _get_correlation_matrix(map, rms_x, rms_y, rms_z, beta)
    if comm == None:
        comm = MPI.COMM_WORLD
    bunch = Bunch(lattice_simulator.get_lattice().get_reference_particle(),
                  num_macro_particles, num_real_particles, comm);
    dist = Random_distribution(seed, comm)
    populate_6d(dist, bunch, numpy.zeros((6,), 'd'), correlation_matrix)
    return bunch

def get_matched_bunch_transverse_parameters(lattice_simulator,
                                            emit_x, emit_y, rms_z, rms_dpop):

    map = linear_one_turn_map(lattice_simulator)
    alpha, beta = get_alpha_beta(map)
    sigma4, r4 = match_transverse_twiss_emittance([emit_x, emit_y], alpha, beta)
    sigma = range(0, 6)
    sigma[0:4] = sigma4
    beta = lattice_simulator.get_lattice().get_reference_particle().get_beta()
    p = lattice_simulator.get_lattice().get_reference_particle().get_momentum()
    sigma[4] = rms_z / beta
    sigma[5] = rms_dpop
    r = [r4[0], r4[1], 0.0]
    covariance_matrix = get_covariances(sigma, r)
    means = numpy.zeros((6,), 'd')
    return means, covariance_matrix

def generate_matched_bunch_transverse(lattice_simulator, emit_x, emit_y,
                           rms_z, dpop, num_real_particles,
                           num_macro_particles, seed=0, comm=None):

    means, covariance_matrix = \
        get_matched_bunch_transverse_parameters(lattice_simulator,
                                                emit_x, emit_y, rms_z, dpop)
    if comm == None:
        comm = MPI.COMM_WORLD
    bunch = Bunch(lattice_simulator.get_lattice().get_reference_particle(),
                  num_macro_particles, num_real_particles, comm)
    dist = Random_distribution(seed, comm)
    populate_6d(dist, bunch, means, covariance_matrix)
    return bunch

def generate_matched_bunch_uniform_longitudinal(lattice_simulator, emit_x,
                           emit_y, length_cdt, dpop, num_real_particles,
                           num_macro_particles, seed=0, comm=None):

    means, covariance_matrix = \
        get_matched_bunch_transverse_parameters(lattice_simulator,
                                                emit_x, emit_y, 1.0, dpop)
    if comm == None:
        comm = MPI.COMM_WORLD
    bunch = Bunch(lattice_simulator.get_lattice().get_reference_particle(),
                  num_macro_particles, num_real_particles, comm)
    dist = Random_distribution(seed, comm)
    populate_transverse_gaussian(dist, bunch, means, covariance_matrix,
                                 length_cdt)
    return bunch
