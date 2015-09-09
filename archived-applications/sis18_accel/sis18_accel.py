#!/usr/bin/env python
import sys
from math import sqrt
import synergia
from sis18_accel_options import opts
from mpi4py import MPI
import numpy as np
import tables
import os

import read_bunch
from ramp_module import Ramp_actions
#================================================
# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        raise RuntimeError, "map is unstable"

    mu = np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)
#================================================


# get_fd_quads(lattice) reads input lattice and returns
# ( [list of focussing quad elements], [list of defocussing quad elements] )
def get_fd_quads(lattice):
    f_quads = []
    d_quads = []
    for elem in lattice.get_elements():
        if elem.get_type() == "quadrupole":
            etype = elem.get_string_attribute("type").lower()
            if (etype == "qd1"):
                f_quads.append(elem)
            elif (etype == "qd2"):
                d_quads.append(elem)
    return (f_quads, d_quads)
#================================================

try:
    logger = synergia.utils.Logger(0)

    # this is the communicator object that will be used for MPI operations
    comm = synergia.utils.Commxx()
    myrank = comm.get_rank()
    mpisize = comm.get_size()
    verbose = opts.verbosity>0

#=================  Setting up the lattice ==================================

    lattice = synergia.lattice.Mad8_reader().get_lattice("machine", "sis18-6-accel.mad")
    lattice_length = lattice.get_length()
    print >>logger, "Lattice length: ", lattice.get_length()

    f_quads, d_quads = get_fd_quads(lattice)
    print >>logger, "There are ", len(f_quads), " focussing quadrupoles"
    print >>logger, "There are ", len(d_quads), " defocussing quadrupoles"

    # set k values for focussing and defocussing quads
    kq = opts.kqf
    for elem in f_quads:
        elem.set_double_attribute("k1", kq)
    kq = opts.kqd
    for elem in d_quads:
        elem.set_double_attribute("k1", kq)

    reference_particle = lattice.get_reference_particle()
    energy = reference_particle.get_total_energy()
    beta = reference_particle.get_beta()
    gamma = reference_particle.get_gamma()

    h = 1.0
    # calculate RF frequency for informational purposes
    freq = h * beta * synergia.foundation.pconstants.c*1.0e-6/lattice_length

    print >>logger, "lattice length: ", lattice_length
    print >>logger, "energy: ", energy
    print >>logger, "beta: ", beta
    print >>logger, "gamma: ", gamma
    print >>logger, "RF frequency: ", freq, " MHz"

    lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                              opts.map_order)
    # adjust tunes before setting sextupole and RF

    compaction_factor = lattice_simulator.get_momentum_compaction()
    slip_factor = lattice_simulator.get_slip_factor()
    print >>logger, "compaction factor (alfa): ", compaction_factor
    print >>logger, "slip factor: ", slip_factor

    if opts.xtune:
        (orig_xtune, orig_ytune) = lattice_simulator.get_both_tunes()
        new_xtune = opts.xtune
        print >>logger, "Adjusting x tune to: ", new_xtune

        lattice_simulator.adjust_tunes(new_xtune, orig_ytune, f_quads, d_quads,
                                       opts.tune_tolerance)

    # set sextupole if requested
    if opts.k2l != 0:

        for elem in lattice.get_elements():
            if elem.get_type() == "multipole" and elem.get_name()=="s":
                print >>logger, "setting k2l to: ", opts.k2l, " in element: ", elem.get_name()
                elem.set_double_attribute("k2l", opts.k2l)

    # Set the same aperture radius for all elements,  RF frequency
    # and set the correct propagator
    if opts.full_maps:
        # use maps for propagation
        print >>logger, "Using maps for propagation"
    else:
        # use CHEF symplectic tracking
        print >>logger, "Using chef for propagation"

    for elem in lattice.get_elements():

        # set RF frequency so the lattice simulator can find the
        # bucket length
        if elem.get_type() == "rfcavity":
            elem.set_double_attribute("freq", freq)

        elem.set_double_attribute("circular_aperture_radius", opts.radius)

        if opts.full_maps:
            elem.set_string_attribute("extractor_type", "chef_map")
        else:
            elem.set_string_attribute("extractor_type", "chef_propagate")

    lattice_simulator.update()

    if opts.xml_save_lattice:
        synergia.lattice.xml_save_lattice(lattice, "sis18-6.xml")

#=================  Finished setting up the lattice ==================================

#================= get lattice properties ==========================

    map = lattice_simulator.get_linear_one_turn_map()
    if verbose:
        print >>logger, "one turn map"
        print >>logger, np.array2string(map, max_line_width=200)

    [l, v] = np.linalg.eig(map)

    if verbose:
        print >>logger, "eigenvalues: "
        for z in l:
            print >>logger, "|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi)

    [ax, bx, qx] = map2twiss(map[0:2,0:2])
    [ay, by, qy] = map2twiss(map[2:4, 2:4])
    [az, b_cdt, qz] = map2twiss(map[4:6,4:6])

    if verbose:
        print >>logger, "Lattice parameters (assuming uncoupled map)"
        print >>logger, "alpha_x: ", ax, " alpha_y: ", ay
        print >>logger, "beta_x: ", bx, " beta_y: ", by
        print >>logger, "q_x: ", qx, " q_y: ", qy, "q_s: ", qz
        print >>logger, "beta_cdt: ", b_cdt, " beta_longitudinal: ", b_cdt*beta

#=================  Finish getting lattice properties  ==================================

#================= Either generate or read in a bunch ===================================

    if opts.generate_bunch:

        # emittance in options in 2 sigma
        emitx = opts.emitx
        emity = opts.emity
        stdx = sqrt(bx*emitx/4.0)
        stdy = sqrt(by*emity/4.0)
        dpop = opts.dpop
        std_cdt = dpop*b_cdt

        if verbose:
            print >>logger, "dpop",dpop
            print >>logger, "std_cdt: ",std_cdt

        if verbose:
            print >>logger, "Beam parameters for bunch generation"
            print >>logger, "invariant emittance_H: ", opts.emitx, " [m-rad]"
            print >>logger, "invariant emittance_V: ", opts.emity, " [m-rad]"
            print >>logger, "Delta-p/p RMS: ", dpop
            print >>logger, "RMS x: ", stdx, " [m]"
            print >>logger, "RMS y: ", stdy, " [m]"
            print >>logger, "RMS z: ", std_cdt*beta, " [m], RMS c*dt: ", std_cdt, " [m], RMS dt: ", std_cdt*1.0e9/synergia.foundation.pconstants.c, " [ns]"
            print >>logger
            if opts.test_particles:
                print >>logger, "Test particles ARE included"
            else:
                print >>logger, "Test particles ARE NOT included"


        if opts.matching == "normalform":

            # get actions that give an equilibrium distribution with desired 2nd moments
            actions = lattice_simulator.get_stationary_actions(stdx, stdy, std_cdt)

            if (actions[0]<0.0) or (actions[1]<0.0) or (actions[2]<0.0):
                raise RuntimeError, "get_stationary_actions can't satisfy requested moments"

            # generate a periodic bunch
            bunch = synergia.bunch.Bunch(
                lattice.get_reference_particle(),
                opts.macro_particles, opts.real_particles, comm,
                lattice_simulator.get_bucket_length(), 0)

            seed = opts.seed
            dist = synergia.foundation.Random_distribution(seed, comm)

            half_bucket_length = lattice_simulator.get_bucket_length()/2.0
            if verbose:
                print >>logger, "half_bucket_length: ", half_bucket_length, " [m]"

            # generation limits are really given in cdt
            synergia.simulation.populate_6d_stationary_clipped_longitudinal_gaussian(
                dist, bunch, actions, -half_bucket_length/beta, half_bucket_length/beta, lattice_simulator)
        elif opts.matching == "6dlinear":
            # match xrms, yrms, dpop
            consigma = opts.cutoffnsigma*np.ones([6])
            consigma = np.array([6.0, 6.0, 6.0, 6.0, opts.cutoffnsigma, opts.cutoffnsigma])
            bunch = synergia.optics.generate_matched_bunch(
                lattice_simulator,
                stdx, stdy, dpop,
                opts.real_particles,
                opts.macro_particles,
                rms_index=[0,2,5],
                periodic=True,
                seed=opts.seed, nsigma=consigma)
        else:
            raise RuntimeError, "Invalid matching procedure %d specified"%opts.matching

        # create 201 test particles with 0 particle and 50 horizontal and vertical particles at .1 sigma
        # to 5 sigma, plus another 50 x and y at exact negative coordintes


        # I have to distribute these particles across all processors

        # I have 101 particles to distribute over mpisize processors.
        # To preserve symmetry, I'll generate 201 particle symmetrically.
        # First I'll make the particles (all processors)
        test_particles = np.zeros((201, 7),'d')
        # 50 horizontal test particles from .1 to 5 sigma
        # particle 0 is OK as is

        curpos = 1
        tstp = opts.test_step
        for tsig in range(1,51):
            test_particles[curpos,0:6] = 0.0
            test_particles[curpos,0] = (tstp*tsig) * opts.x_test
            # now we start with particles on the x axis
            #test_particles[curpos,1] = -ax/bx * ((tstp*tsig) * opts.x_test)
            # possibly have particles start with z offset
            test_particles[curpos,4] = opts.z_offset/beta
            test_particles[curpos,6] = float(curpos)
            curpos = curpos + 1
            
        # 50 vertical test particles from .1 to 5 sigma
        for tsig in range(1,51):
            test_particles[curpos,0:6] = 0.0
            test_particles[curpos,2] = (tstp*tsig) * opts.y_test
            # now we start with particles on the y axis
            #test_particles[curpos,3] = -ay/by * ((tstp*tsig) * opts.y_test)
            # possibly have particles start with z offset
            test_particles[curpos,4] = opts.z_offset/beta
            test_particles[curpos,6] = float(curpos)
            curpos = curpos + 1

        # Now the symmetric particles on the opposite side
        for tsig in range(1,51):
            test_particles[curpos,0:6] = 0.0
            test_particles[curpos,0] = -(tstp*tsig) * opts.x_test
            #test_particles[curpos,1] = ax/bx * ((tstp*tsig) * opts.x_test)
            # possibly have particles start with z offset
            test_particles[curpos,4] = opts.z_offset/beta
            test_particles[curpos,6] = float(curpos)
            curpos = curpos + 1
            
        # 50 vertical test particles from .1 to 5 sigma
        for tsig in range(1,51):
            test_particles[curpos,0:6] = 0.0
            test_particles[curpos,2] = -(tstp*tsig) * opts.y_test
            #test_particles[curpos,3] = ay/by * ((tstp*tsig) * opts.y_test)
            # possibly have particles start with z offset
            test_particles[curpos,4] = opts.z_offset/beta
            test_particles[curpos,6] = float(curpos)
            curpos = curpos + 1


        if curpos != 201:
            raise RuntimeError, "wrong number of test particles"

        # distribute the test particles to appropriate processors
        # each processor will have a few test particles at the beginning
        # of its local_particles.  First figure out how many are
        # on each processor

        ntest_particles = curpos
        offsets,counts = synergia.utils.decompose_1d_raw(mpisize, ntest_particles)

        # move the correct test particles into the local array
        local_particles = bunch.get_local_particles()

        if opts.test_particles:
            # now it gets tricky.  I'm going to move the centroid of the
            # bunch to put it at 0, but I need to calculate the centroid taking
            # into account the test particles, so first 0 the test particle
            # locations before calculating centroid.  The test particle centroid is
            # 0 by construction
        
            for ppart in range(counts[myrank]):
                local_particles[ppart, 0:6] = 0.0
                local_particles[ppart, 6] = test_particles[offsets[myrank], 6]

            # now calculate and subtract centroid.  The test particles I am going
            # to be adding later have mean 0, but the mean I subtract has to be
            # weighted by the total number of particles-201
            total_num = bunch.get_total_num()
            reweight = total_num/float(total_num - ntest_particles)
            means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)*reweight
            for k in range(6):
                local_particles[:,k] = local_particles[:,k] - means[k]
        
            # Now put in my test particles
            for ppart in range(counts[myrank]):
                local_particles[ppart,:] = test_particles[offsets[myrank]+ppart,:]

        # not using test particles
        else:
            means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
            for k in range(6):
                local_particles[:,k] = local_particles[:,k] - means[k]
            

    else:
        # read the bunch from a file of MAD-X particle data
        particles_file = opts.particles_file
        print >>logger, "Reading bunch particles from file: ", particles_file

        bunch = read_bunch.read_bunch(opts.particles_file, lattice.get_reference_particle(), lattice_simulator.get_bucket_length(), comm, verbose)


    if verbose:
        read_bunch.print_bunch_stats(bunch)

#============================ Finished acquiring a bunch =====================

#================== set up the output diagnostics by creating the bunch_simulator =============

    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)


    if opts.accel:
        step_diag_prefix = "step_accel_"
        turn_diag_prefix = "turn_accel_"
    else:
        step_diag_prefix = "step_noaccel_"
        turn_diag_prefix  = "turn_noaccel_"

    if opts.step_full2:
        if opts.scratch:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2(step_diag_prefix+"full2.h5", opts.scratch))
        else:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2(step_diag_prefix+"step_full2.h5"))
        print >>logger, "saving full2 diagnostics each step"

    if opts.turn_full2:
        if opts.scratch:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2(turn_diag_prefix+"full2.h5", opts.scratch))
        else:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2(turn_diag_prefix+"full2.h5"))
        print >>logger, "saving full2 diagnostics each turn"

    if opts.step_tracks > 0:
        if opts.scratch:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_bulk_track(step_diag_prefix+"tracks.h5", opts.step_tracks, 0, opts.scratch))
        else:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_bulk_track(step_diag_prefix+"tracks.h5", opts.step_tracks, 0))
        print >>logger, "saving particle track data for ", opts.step_tracks, " particls"

    if opts.turn_tracks > 0:
        if opts.scratch:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track(turn_diag_prefix+"tracks.h5", opts.turn_tracks, 0, opts.scratch))
        else:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track(turn_diag_prefix+"tracks.h5", opts.turn_tracks, 0))
        print >>logger, "saving particle track data for ", opts.turn_tracks, " particls"

    if opts.step_particles:
        if opts.scratch:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_particles(step_diag_prefix+"particles.h5", 0, 0, opts.scratch))
        else:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_particles(step_diag_prefix+"particles.h5", 0, 0))
        print >>logger, "saving turn-by-turn particle data each step"

    if opts.turn_particles:
        if opts.scratch:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles(turn_diag_prefix+"particles.h5",0, 0, opts.scratch), opts.turn_particles_period)
        else:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles(turn_diag_prefix+"particles.h5",0, 0), opts.turn_particles_period)
        print >>logger, "saving turn-by-turn particle data each turn every %d turns"%opts.turn_particles_period

        bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_reference_particle("refpart.h5", ""))

#=========== finished setting up diagnostics =============================

#============================ Turn on acceleration =====================

    for elem in lattice.get_elements():
        if elem.get_type() == "rfcavity":
            if opts.accel:
                print >>logger, "Setting RF cavity for acceleration"
                # for half the maximum voltage, use 30 degrees which is
                # 1/12 of 2 pi (the units of lag) which gives 1 MeV/turn.
                #elem.set_double_attribute("lag", 1.0/12.0)
                # a smaller emergy gain
                elem.set_double_attribute("lag", np.arcsin(0.5*0.1)/(2.0*np.pi))
            else:
                print >>logger, "Setting RF cavity for no acceleration"


    stepper = synergia.simulation.Independent_stepper(lattice, opts.map_order, opts.steps)

    propagator = synergia.simulation.Propagator(stepper)
    propagator.set_concurrent_io(opts.concurrent_io)
    propagator.set_checkpoint_period(opts.checkpointperiod)

    ramp_actions = Ramp_actions(bunch.get_reference_particle().get_total_energy())

    propagator.propagate(bunch_simulator, ramp_actions, opts.turns, opts.maxturns, opts.verbosity)

    sys.exit(0)

except Exception, e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
