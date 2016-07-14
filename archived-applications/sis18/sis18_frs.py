#!/usr/bin/env python
import sys
from math import sqrt
import synergia
from sis18_frs_options import opts
from mpi4py import MPI
import numpy as np
import tables
import os

if opts.dump_efield:
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
            etype = elem.get_string_attribute("type")
            if (etype == "qd1"):
                f_quads.append(elem)
            elif (etype == "qd2"):
                d_quads.append(elem)
    return (f_quads, d_quads)
#================================================


# load the particles that will be used for the simulation
# The particles file is a text file with particle coordinates
# defined with the MAD-X conventions: X PX Y PY T PT
# Read this in using numpy's loadtxt command
# particle coordinates are converted to Synergia conventions

# input arguments:
#    particles_file: the file name
#    reference particle: the lattice reference particle for kinematic conversions
#    comm: the Commxx communicator object for this bunch

#  
def read_bunch(particles_file, refpart, bucket_length, comm, verbose=False):
    name,extension = os.path.splitext(particles_file)
    if extension == ".txt":
        return read_txt_particles(particles_file, refpart, bucket_length, comm, verbose)
    elif extension == ".h5":
        return read_h5_particles(particles_file, refpart, bucket_length, comm, verbose)

#====================================================================

def read_txt_particles(particles_file, refpart, bucket_length, comm, verbose):
    four_momentum = refpart.get_four_momentum()
    pmass = four_momentum.get_mass()
    E_0 = four_momentum.get_total_energy()
    p0c = four_momentum.get_momentum()

    myrank = comm.get_rank()
    mpisize = comm.get_size()
    
    if myrank==0:
        print "Loading particles from txt file: ", particles_file

    if myrank == 0:
        particles = np.loadtxt(particles_file)
        num_total_particles = particles.shape[0]
        # broadcast num particles to all nodes
        MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    else:
        num_total_particles = None
        num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    if myrank == 0:
        # make sure the data has the correct shape, either [n,6] without
        # particles IDs or [n,7] with particle IDs.
        if (particles.shape[1] != 6) and (particles.shape[1] != 7):
            raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)
        
        # numpy manipulations to convert kinematics
        # convert MAD-X T=-c*dt to Synergia c*ct
        particles[:,4] = -particles[:,4]
        # convert MAD-X Delta-E/pc to Synergia delta-p/p
        # sqrt(((dE/p0c)+(E0/p0c))**2 - (m/p0c)**2) - (p0c/p0c)
        m_over_pc = pmass/p0c
        E_0_over_pc = E_0/p0c
        particles[:,5] = np.sqrt( (particles[:,5] + E_0_over_pc) *
                                  (particles[:,5] + E_0_over_pc) - m_over_pc**2 ) - 1.0
        

        # if there are no IDs, append particle ID column
        if particles.shape[1] != 7:
            particles_w_id = np.column_stack((particles,
                                              np.arange(num_total_particles, dtype='d')))
        else:
            particles_w_id = particles
            
            if myrank == 0:
                print "Read ", num_total_particles, " particles"


    # create a bunch with the correct number of macro particles
    bunch = synergia.bunch.Bunch(
        refpart,
        num_total_particles, opts.real_particles, comm,
        bucket_length)

    local_num = bunch.get_local_num()
    local_particles = bunch.get_local_particles()

    # Each processor will have a possibly different number of local particles.
    # rank 0 has to find out how many each of them has and distribute them
    n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    if myrank == 0:
        # copy in my particles
        this_rank_start = 0
        local_particles[:,:] = particles_w_id[0:local_num, :]
        this_rank_start += local_num
        # send particles out to other ranks
        for r in range(1, mpisize):
            this_rank_end = this_rank_start+n_particles_by_proc[r]
            MPI.COMM_WORLD.send(obj=particles_w_id[this_rank_start:this_rank_end, :],
                                dest=r)
            this_rank_start += n_particles_by_proc[r]
    else:
        # I'm not rank 0.  Receive my particles
        lp = MPI.COMM_WORLD.recv(source=0)
        local_particles[:,:] = lp[:,:]
    return bunch

#==========================================================

def read_h5_particles(particles_file, refpart, bucket_length, comm, verbose):
    four_momentum = refpart.get_four_momentum()
    pmass = four_momentum.get_mass()
    E_0 = four_momentum.get_total_energy()
    p0c = four_momentum.get_momentum()

    myrank = comm.get_rank()
    mpisize = comm.get_size()
    
    if myrank==0:
        print "Loading particles from h5 file: ", particles_file

    if myrank == 0:
        h5 = tables.openFile(particles_file)
        # use explicit int conversion otherwise there seems to
        # be a typepython->C++ type  mismatch of numpy.int64->int
        num_total_particles = int(h5.root.particles.shape[0])
        if verbose:
            print "Total of  ", num_total_particles, " particles from file"
        # broadcast num particles to all nodes
        MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    else:
        num_total_particles = None
        num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    if myrank == 0:
        particles = h5.root.particles
        # make sure the data has the correct shape, either [n,6] without
        # particles IDs or [n,7] with particle IDs.
        if (particles.shape[1] != 7):
            raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)


    # create a bunch with the correct number of macro particles
    bunch = synergia.bunch.Bunch(
        refpart,
        num_total_particles, opts.real_particles, comm,
        bucket_length, 0)

    local_num = bunch.get_local_num()
    local_particles = bunch.get_local_particles()

    # Each processor will have a possibly different number of local particles.
    # rank 0 has to find out how many each of them has and distribute them
    n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    if myrank == 0:
        # copy in my particles
        this_rank_start = 0
        local_particles[:,:] = particles[0:local_num, :]
        this_rank_start += local_num
        # send particles out to other ranks
        for r in range(1, mpisize):
            this_rank_end = this_rank_start+n_particles_by_proc[r]
            MPI.COMM_WORLD.send(obj=particles[this_rank_start:this_rank_end, :],
                                dest=r)
            this_rank_start += n_particles_by_proc[r]
    else:
        # I'm not rank 0.  Receive my particles
        lp = MPI.COMM_WORLD.recv(source=0)
        local_particles[:,:] = lp[:,:]

    return bunch

#================================================================

def print_bunch_stats(bunch):
    coord_names = ("x", "xp", "y", "yp", "c*dt", "dp/p")

    means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
    stds = synergia.bunch.Core_diagnostics().calculate_std(bunch, means)
    if myrank == 0:
        print "%20s   %20s   %20s"%("coord","mean","rms")
        print "%20s   %20s   %20s"%("====================",
                                    "====================",
                                    "====================")
        for i in range(6):
            print "%20s   %20.12e   %20.12e"%(coord_names[i], means[i], stds[i])

#=========================================================
try:
    # this is the communicator object that will be used for MPI operations
    comm = synergia.utils.Commxx()
    myrank = comm.get_rank()
    mpisize = comm.get_size()
    verbose = opts.verbosity>0

#=================  Setting up the lattice ==================================

    lattice = synergia.lattice.Mad8_reader().get_lattice("machine", "sis18-6.mad")
    lattice_length = lattice.get_length()

    f_quads, d_quads = get_fd_quads(lattice)
    if myrank == 0:
        print "There are ", len(f_quads), " focussing quadrupoles"
        print "There are ", len(d_quads), " defocussing quadrupoles"

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
    # recalculate RF frequency so it is numerically consistent with
    # lattice length within machine precision
    freq = h * beta * synergia.foundation.pconstants.c*1.0e-6/lattice_length

    if myrank == 0:
        print "lattice length: ", lattice_length
        print "energy: ", energy
        print "beta: ", beta
        print "gamma: ", gamma
        print "RF frequency: ", freq, " MHz"

    lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                              opts.map_order)
    # adjust tunes before setting sextupole and RF

    compaction_factor = lattice_simulator.get_momentum_compaction()
    slip_factor = lattice_simulator.get_slip_factor()
    print "compaction factor (alfa): ", compaction_factor
    print "slip factor: ", slip_factor

    if opts.xtune:
        (orig_xtune, orig_ytune) = lattice_simulator.get_both_tunes()
        new_xtune = opts.xtune
        if myrank == 0:
            print "Adjusting x tune to: ", new_xtune

        lattice_simulator.adjust_tunes(new_xtune, orig_ytune, f_quads, d_quads,
                                       opts.tune_tolerance)

    # set sextupole if requested
    if opts.k2l != 0:

        for elem in lattice.get_elements():
            if elem.get_type() == "multipole" and elem.get_name()=="s":
                if myrank == 0:
                    print "setting k2l to: ", opts.k2l, " in element: ", elem.get_name()
                elem.set_double_attribute("k2l", opts.k2l)

    # Set the same aperture radius for all elements and fix RF frequency
    for elem in lattice.get_elements():

        # set RF frequency so the lattice simulator can find the
        # bucket length
        if elem.get_type() == "rfcavity":
            elem.set_double_attribute("freq", freq)
            #elem.set_double_attribute("lag", 0.0)
            #elem.set_double_attribute("volt", 0.0)

        elem.set_double_attribute("circular_aperture_radius", opts.radius)

        # set the propagation method.
        # extractor_type may be "chef_maps", "chef_mixed", or
        #    chef_propagate
        if opts.use_maps == "all":
            elem.set_string_attribute("extractor_type", "chef_map")
        elif opts.use_maps == "none":
            elem.set_string_attribute("extractor_type", "chef_propagate")
        elif opts.use_maps == "nonrf":
            elem.set_string_attribute("extractor_type", "chef_mixed")
        elif opts.use_maps == "onlyrf":
            if elem.get_type() == "rfcavity":
                elem.set_string_attribute("extractor_type", "chef_map")
            else:
                elem.set_string_attribute("extractor_type", "chef_propagate")
        elif opts.use_maps == "libff":
            elem.set_string_attribute("extractor_type", "libff")
        else:
            raise RuntimeError, "bad options for use_maps: %d"%opts.use_maps

    lattice_simulator.update()

    if opts.xml_save_lattice:
        synergia.lattice.xml_save_lattice(lattice, "sis18-6.xml")

#=================  Finished setting up the lattice ==================================

#================= get lattice properties ==========================

    map = lattice_simulator.get_linear_one_turn_map()
    if verbose and (myrank==0):
        print "one turn map from synergia2.5 infrastructure"
        print np.array2string(map, max_line_width=200)

    [l, v] = np.linalg.eig(map)

    if verbose and myrank==0:
        print "eigenvalues: "
        for z in l:
            print "|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi)

    [ax, bx, qx] = map2twiss(map[0:2,0:2])
    [ay, by, qy] = map2twiss(map[2:4, 2:4])
    [az, b_cdt, qz] = map2twiss(map[4:6,4:6])

    if verbose and (myrank == 0):
        print "Lattice parameters (assuming uncoupled map)"
        print "alpha_x: ", ax, " alpha_y: ", ay
        print "beta_x: ", bx, " beta_y: ", by
        print "q_x: ", qx, " q_y: ", qy, "q_s: ", qz
        print "beta_cdt: ", b_cdt, " beta_longitudinal: ", b_cdt*beta

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

        if verbose and (myrank == 0):
            print "dpop",dpop
            print "std_cdt: ",std_cdt



        if verbose and (myrank == 0):
            print "Beam parameters for bunch generation"
            print "invariant emittance_H: ", opts.emitx, " [m-rad]"
            print "invariant emittance_V: ", opts.emity, " [m-rad]"
            print "Delta-p/p RMS: ", dpop
            print "RMS x: ", stdx, " [m]"
            print "RMS y: ", stdy, " [m]"
            print "RMS z: ", std_cdt*beta, " [m], RMS c*dt: ", std_cdt, " [m], RMS dt: ", std_cdt*1.0e9/synergia.foundation.pconstants.c, " [ns]"
            print
            if opts.test_particles:
                print "Test particles ARE included"
            else:
                print "Test particles ARE NOT included"


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
            if verbose and (myrank == 0):
                print "half_bucket_length: ", half_bucket_length, " [m]"

            # generation limits are really given in cdt
            synergia.simulation.populate_6d_stationary_clipped_longitudinal_gaussian(
                dist, bunch, actions, -half_bucket_length/beta, half_bucket_length/beta, lattice_simulator)
        elif opts.matching == "6dlinear":
            # match xrms, yrms, dpop
            consigma = opts.cutoffnsigma*np.ones([6])
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

        DEBUG_DISTRIBUTION = False

        if DEBUG_DISTRIBUTION and myrank==0:
            print "offsets: ", offsets
            print "counts: ", counts

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

            if DEBUG_DISTRIBUTION and myrank == 0:
                print "statistics before subtracting centroid"
            if DEBUG_DISTRIBUTION:
                print_bunch_stats(bunch)

            # now calculate and subtract centroid.  The test particles I am going
            # to be adding later have mean 0, but the mean I subtract has to be
            # weighted by the total number of particles-201
            total_num = bunch.get_total_num()
            reweight = total_num/float(total_num - ntest_particles)
            means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)*reweight
            for k in range(6):
                local_particles[:,k] = local_particles[:,k] - means[k]
        
            if DEBUG_DISTRIBUTION and myrank == 0:
                print "statistics after subtracting centroid"
            if DEBUG_DISTRIBUTION:
                print_bunch_stats(bunch)
        
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
        if myrank == 0:
            print "Reading bunch particles from file: ", particles_file

        bunch = read_bunch(opts.particles_file, lattice.get_reference_particle(), lattice_simulator.get_bucket_length(), comm, verbose)


    if verbose:
        print_bunch_stats(bunch)

#============================ Finished acquiring a bunch =====================

#============================ Set up the stepper and space charge solver =====

    requested_stepper = opts.stepper
    if myrank == 0:
        print "requested_stepper: ",  requested_stepper
        print "use maps for: ", opts.use_maps
#====================== make space charge solver if needed ==================

    if opts.spacecharge:
        solver = opts.solver

        # space charge only works with the split operator stepper, or soelements 
        if (requested_stepper != "splitoperator") and (requested_stepper != "soelements"):
            requested_stepper = "soelements"
            if myrank == 0:
                print "requested stepper changed to soelements for space charge"

        gridx = opts.gridx
        gridy = opts.gridy
        gridz = opts.gridz
        grid = [gridx, gridy, gridz]

        if myrank == 0:
            print "grid: ", grid

        if opts.comm_divide:
            sc_comm = synergia.utils.Commxx_divider(opts.comm_divide, False)
        else:
            sc_comm = synergia.utils.Commxx(True)

        #sc_comm = synergia.utils.Commxx(True)
        if solver == "2dopen-hockney":
            coll_operator = synergia.collective.Space_charge_2d_open_hockney(sc_comm, grid)
        elif solver == "3dopen-hockney":
            # full signature for 3d_open_hockney constructor is
            # comm, grid, long_kicks, z_periodic, period, grid_entire_period,
            # nsigma

            coll_operator = synergia.collective.Space_charge_3d_open_hockney(sc_comm, grid, opts.long_kicks, False, 0.0, False, opts.nsigma)
        elif solver == "2dbassetti-erskine":
            coll_operator = synergia.collective.Space_charge_2d_bassetti_erskine()
        else:
            raise RuntimeError, "requested space charge operator %s invalid.  Must be either 2dopen-hockney or 3dopen-hockney"%opts.solver

        if myrank == 0:
            print "Using space charge solver ", solver
            print "Grid: ", gridx, " x ", gridy, " x ", gridz

    else:
        coll_operator = synergia.simulation.Dummy_collective_operator("stub")
        if myrank == 0:
            print "No space charge solver used"

#==================== Finished space charge solver ======================

#==================== set up the stepper ===============================

    if requested_stepper == "splitoperator":

        if myrank == 0:
            print "Using split-operator stepper with ", opts.steps, " steps/turn"

        stepper = synergia.simulation.Split_operator_stepper(
            lattice, opts.map_order, coll_operator, opts.steps)

    elif requested_stepper == "soelements":

        if myrank == 0:
            print "Using split-operator stepper elements with ", opts.steps, " steps/element"

        stepper = synergia.simulation.Split_operator_stepper_elements(
            lattice, opts.map_order, coll_operator, opts.steps)

    elif requested_stepper == "independent":

        if myrank == 0:
            print "Using independent-operator stepper with ", opts.steps, " steps/turn"

        stepper = synergia.simulation.Independent_stepper(
            lattice, opts.map_order, opts.steps)

    elif requested_stepper == "elements":

        if myrank == 0:
            print "Using step-by-elements-operator stepper with ", opts.steps, " steps/element"

        stepper = synergia.simulation.Independent_stepper_elements(
            lattice, opts.map_order, opts.steps)

    else:
        raise RuntimeError, "stepper %s invalid, must be either 'splitoperator', 'independent' or 'elements'"%requested_stepper

#================== finished setting up the stepper ========================

#================== set up the output diagnostics by creating the bunch_simulator =============

    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)


    if opts.step_full2:
        if opts.scratch:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2("step_full2.h5", opts.scratch))
        else:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2("step_full2.h5"))
        if myrank == 0:
            print "saving full2 diagnostics each step"
    if opts.turn_full2:
        if opts.scratch:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("turn_full2.h5", opts.scratch))
        else:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("turn_full2.h5"))
        if myrank == 0:
            print "saving full2 diagnostics each turn"

    if opts.step_tracks > 0:
        if opts.scratch:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_bulk_track("step_tracks.h5", opts.step_tracks, 0, opts.scratch))
        else:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_bulk_track("step_tracks.h5", opts.step_tracks, 0))
        if myrank == 0:
            print "saving particle track data for ", opts.step_tracks, " particls"

    if opts.turn_tracks > 0:
        if opts.scratch:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("turn_tracks.h5", opts.turn_tracks, 0, opts.scratch))
        else:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("turn_tracks.h5", opts.turn_tracks, 0))
        if myrank == 0:
            print "saving particle track data for ", opts.turn_tracks, " particls"

    if opts.step_particles:
        if opts.scratch:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_particles("step_particles.h5", 0, 0, opts.scratch))
        else:
            bunch_simulator.add_per_step(synergia.bunch.Diagnostics_particles("step_particles.h5", 0, 0))
        if myrank == 0:
            print "saving turn-by-turn particle data each step"

    if opts.turn_particles:
        if opts.scratch:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("turn_particles.h5",0, 0, opts.scratch), opts.turn_particles_period)
        else:
            bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("turn_particles.h5",0, 0), opts.turn_particles_period)
        if myrank == 0:
            print "saving turn-by-turn particle data each turn every %d turns"%opts.turn_particles_period


    if opts.per_operator_full2:
        bunch_simulator.get_diagnostics_actions().add_per_operator(synergia.bunch.Diagnostics_full2("operator_full2.h5"))

    if opts.per_operator_tracks:
        bunch_simulator.get_diagnostics_actions().add_per_operator(synergia.bunch.Diagnostics_bulk_track("operator_tracks.h5", opts.per_operator_tracks, 0))

    if opts.per_operator_particles:
        bunch_simulator.get_diagnostics_actions().add_per_operator(synergia.bunch.Diagnostics_particles("operator_particles.h5"))

#=========== finished setting up diagnostics ===========================================================

#=========== put it all together in the propagator and go!

    # save initial  particle distributions no matter what
    sis18_init = synergia.bunch.Diagnostics_particles("sis18_initial.h5")
    sis18_init.set_bunch(bunch)
    sis18_init.update_and_write()

    if opts.dump_efield:
        ramp_actions = Ramp_actions(coll_operator)

    propagator = synergia.simulation.Propagator(stepper)
    propagator.set_checkpoint_period(opts.checkpointperiod)
    propagator.set_concurrent_io(opts.concurrent_io)

    if opts.dump_efield:
        propagator.propagate(bunch_simulator, ramp_actions, opts.turns, opts.maxturns, opts.verbosity)
    else:
        propagator.propagate(bunch_simulator, opts.turns, opts.maxturns, opts.verbosity)

except Exception, e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
