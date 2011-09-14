#!/usr/bin/env python
import sys
import numpy

from one_turn_map import linear_one_turn_map
from mpi4py import MPI
from synergia.bunch import Bunch, populate_6d, populate_transverse_gaussian, populate_transverse_KV_GaussLong, populate_two_particles
from synergia.foundation import Random_distribution, pconstants
from synergia.optics.one_turn_map import linear_one_turn_map
from math import acos, pi, sin, sqrt


def _get_correlation_matrix(linear_map,arms,brms,crms,beta,rms_index=[0,2,4],print_emittances=True):
    '''here are 3 rms input parameters,arms, brms, crms, which corresponds to  indices rms _index[0], rms _index[1], rms _index[2] 
        example: rms_index=[0,2,4]==> arms=xrms, brms=yrms, crms=zrms
         units of rms should be  [xrms]=m, [pxrms]=Gev/c, [zrms]=m, [pzrms] = Gev/c,  ''' 
    evals, evect_matrix = numpy.linalg.eig(linear_map)
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
            if numpy.max(abs(sum.imag)) < best:
                best = numpy.max(abs(sum.imag))
                conj = item
        if conj == -1:
            raise RuntimeError, "failed to find a conjugate pair in _get_correlation_matrix"
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
    S=numpy.zeros((3,3),'d')
    for i in range(0,3):
        for j in range(0,3):
            S[i,j]=F[j][rms_index[i],rms_index[i]]
        
    Sinv=numpy.linalg.inv(S)   
    
    
    
    gamma=1./numpy.sqrt(1.-beta*beta)
    pz = gamma * beta *pconstants.mp
    energy=pconstants.mp * gamma
    Cxy=1.
    Cxpyp=1./pz
    Cz=1./beta
    Czp=1./pz
    units=[Cxy,Cxpyp,Cxy,Cxpyp,Cz, Czp] # transform from input units to Chef units
 
    C = numpy.zeros([6, 6], 'd') 
    cd1=arms*units[rms_index[0]]*arms*units[rms_index[0]]
    cd2=brms*units[rms_index[1]]*brms*units[rms_index[1]]   
    cd3=crms*units[rms_index[2]]*crms*units[rms_index[2]]  
    #cd1 = rms_x ** 2
    #cd2 = rms_y ** 2
    #rms_cdt = rms_z / beta   
    #cd3 = rms_cdt ** 2
    
    for i in range(0, 3):
        C += F[i] * (Sinv[i, 0] * cd1 + Sinv[i, 1] * cd2 + Sinv[i, 2] * cd3)
        
          
    if print_emittances:
       
        #emitx=sqrt(C[0,0]*C[1,1]-C[0,1]*C[1,0])*pz
        #emity=sqrt(C[2,2]*C[3,3]-C[2,3]*C[3,2])*pz
        #emitz=sqrt(C[4,4]*C[5,5]-C[4,5]*C[5,4])*pz*beta
        emitx=sqrt(C[0,0]*C[1,1]-C[0,1]*C[1,0])/units[0]/units[1]
        emity=sqrt(C[2,2]*C[3,3]-C[2,3]*C[3,2])/units[2]/units[3]
        emitz=sqrt(C[4,4]*C[5,5]-C[4,5]*C[5,4])/units[4]/units[5]
        if MPI.COMM_WORLD.Get_rank() ==0: 
            print "************ BEAM MATCHED PARAMETERS *****************"
            print "*    emitx=", emitx, " meters*GeV/c   =", emitx/pz, " meters*rad =", emitx/pz/pi, " pi*meters*rad"
            print "*    emity=", emity, " meters*GeV/c   =", emity/pz, " meters*rad =", emity/pz/pi, " pi*meters*rad"
            print "*    emitz=", emitz, " meters*GeV/c =", emitz*1.e9/(pconstants.c), " eV*s =" , emitz*beta*beta*energy/pz, " meters*GeV"
            print " "  
            print "*    90%emitx=",  4.605*pi*emitx/pz,"  meters*rad =", 4.605*emitx/pz, " pi*meters*rad"
            print "*    90%emity=",  4.605*pi*emity/pz, " meters*rad =", 4.605*emity/pz, " pi*meters*rad"
            print "*    90%emitz=",  4.605*pi*emitz*1.e9/(pconstants.c), " eV*s"
            print " "
            print " "  
            print "*    95%emitx=",  5.991*pi*emitx/pz,"  meters*rad =", 5.991*emitx/pz, " pi*meters*rad"
            print "*    95%emity=",  5.991*pi*emity/pz, " meters*rad =", 5.991*emity/pz, " pi*meters*rad"
            print "*    95%emitz=",  5.991*pi*emitz*1.e9/(pconstants.c), " eV*s"
            print " "
            print "*    Normalized emitx=",  emitx*gamma*beta/pz, " meters*rad =", emitx*gamma*beta/pz/pi, " pi*meters*rad"
            print "*    Normalized emity=",  emity*gamma*beta/pz, " meters*rad =", emity*gamma*beta/pz/pi, " pi*meters*rad"
            print " " 
            print "*    Normalized 90%emitx=",  4.605*pi*emitx*gamma*beta/pz,"  meters*rad =", 4.605*emitx*gamma*beta/pz, " pi*meters*rad"
            print "*    Normalized 90%emity=",  4.605*pi*emity*gamma*beta/pz, " meters*rad =", 4.605*emity*gamma*beta/pz, " pi*meters*rad"
            print " "  
            print "*    Normalized 95%emitx=",  5.991*pi*emitx*gamma*beta/pz,"  meters*rad =", 5.991*emitx*gamma*beta/pz, " pi*meters*rad"
            print "*    Normalized 95%emity=",  5.991*pi*emity*gamma*beta/pz, " meters*rad =", 5.991*emity*gamma*beta/pz, " pi*meters*rad"
            print " "
            print "*    xrms=",sqrt(C[0,0])/units[0] , " meters"
            print "*    yrms=",sqrt(C[2,2])/units[2], " meters"
            print "*    zrms=",sqrt(C[4,4])/units[4] , " meters=", 1e9*sqrt(C[4,4])/units[4]/pconstants.c/beta," ns  "
           # ,=2.*pi*sqrt(C[4,4])/units[4]/beam_parameters.get_z_length(), " rad 
            print "*    pxrms=",sqrt(C[1,1])/units[1] , " GeV/c,    dpx/p=",sqrt(C[1,1])
            print "*    pyrms=",sqrt(C[3,3])/units[3] , " GeV/c,    dpy/p=",sqrt(C[3,3])
            print "*    pzrms=",sqrt(C[5,5])/units[5] , " GeV/c,    dpz/p=",sqrt(C[5,5])
            print "*    Erms=",sqrt(C[5,5])*beta*beta*energy, " GeV,  deoe=",sqrt(C[5,5])*beta*beta
            #print "" 
            #print "*    bucket length=",beam_parameters.get_z_length(),  " meters =",1e9*beam_parameters.get_z_length()/synergia.PH_MKS_c/beta," ns"
            print "*    pz=",pz, "  GeV/c"
            print "*    total energy=",energy,"GeV,   kinetic energy=", energy-pconstants.mp,"GeV"
            print "****************************************************"
             
            print "One Turn correlation matrix: "
            print numpy.array2string(C, precision=3)

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

def generate_matched_bunch(lattice_simulator, arms,brms,crms,
                           num_real_particles, num_macro_particles, rms_index=[0,2,4],seed=0,
                           comm=None):

    map = linear_one_turn_map(lattice_simulator)
    beta = lattice_simulator.get_lattice().get_reference_particle().get_beta()
    correlation_matrix = _get_correlation_matrix(map, arms,brms,crms,beta, rms_index)
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

def generate_matchedKV_bunch_transverse(lattice_simulator, emitMax,
                           cdt, dpop, num_real_particles,
                           num_macro_particles, seed=0, comm=None):

    map = linear_one_turn_map(lattice_simulator)
    [[ax, ay], [bx, by]] = get_alpha_beta(map)

    if comm == None:
        comm = MPI.COMM_WORLD
    bunch = Bunch(lattice_simulator.get_lattice().get_reference_particle(),
                  num_macro_particles, num_real_particles, comm)
    dist = Random_distribution(seed, comm)
    populate_transverse_KV_GaussLong(dist, bunch, emitMax,
        ax, bx, ay, by, cdt, dpop) 
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

def generate_two_particles(lattice_simulator, coordP1, coordP2, num_real_particles, comm=None):
     if comm == None:
         comm = MPI.COMM_WORLD
     bunch = Bunch(lattice_simulator.get_lattice().get_reference_particle(),
                  2, num_real_particles, comm)
     populate_two_particles(bunch, 
                            coordP1[0], coordP1[1], coordP1[2], coordP1[3], coordP1[4], coordP1[5],
                            coordP2[0], coordP2[1], coordP2[2], coordP2[3], coordP2[4], coordP2[5])
     return bunch

		    
