#!/usr/bin/env python
import sys
import os
import numpy as np
import synergia
import synergia.simulation as SIM
from foborodobo32_options import opts
ET = synergia.lattice.element_type
MT = synergia.lattice.marker_type
#####################################

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        raise RuntimeError("map is unstable")

    mu =np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

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



# mark focussing and defocussing quadrupoles in lattice for adjust tune
def mark_fd_quads(lattice):
    nh = 0
    nv = 0
    for elem in lattice.get_elements():
        if elem.get_type() == ET.quadrupole:
            if elem.get_double_attribute("k1") > 0.0:
                elem.set_marker(MT.h_tunes_corrector)
                elem.reset_marker(MT.v_tunes_corrector)
                nh = nh+1
            elif elem.get_double_attribute("k1") < 0.0:
                elem.reset_marker(MT.h_tunes_corrector)
                elem.set_marker(MT.v_tunes_corrector)
                nv = nv+1
    return (nh, nv)

################################################################################

# mark focussing and defocussing quadrupoles in lattice for adjust tune
def mark_fd_sext(lattice):
    nh = 0
    nv = 0
    for elem in lattice.get_elements():
        if elem.get_type() == ET.sextupole:
            if elem.get_name() == "sf":
                elem.set_marker(MT.h_chrom_corrector)
                elem.reset_marker(MT.v_chrom_corrector)
                nh = nh+1
            elif elem.get_name() == "sd":
                elem.reset_marker(MT.h_chrom_corrector)
                elem.set_marker(MT.v_chrom_corrector)
                nv = nv+1
    return (nh, nv)

################################################################################

# mark focussing and defocussing quadrupoles in lattice for adjust tune
def print_fd_sext(lattice, f):
    nh = 0
    nv = 0
    print('k2 values for horizontal chromaticity correctors')
    for elem in lattice.get_elements():
        if elem.has_marker(MT.h_chrom_corrector):
            if nh != 0:
                print(', ', end='', file=f)
            print(elem.get_double_attribute('k2'), end='', file=f)
            nh = nh + 1
    print(file=f)
    print('k2 values for vertical chromaticity correctors', file=f)
    for elem in lattice.get_elements():
        if elem.has_marker(MT.v_chrom_corrector):
            if nv != 0:
                print(', ', end='', file=f)
            print(elem.get_double_attribute('k2'), end='', file=f)
            nv = nv + 1
    print(file=f)

################################################################################
def print_bunch_stats(bunch, fo):
    coord_names = ("x", "xp", "y", "yp", "c*dt", "dp/p")

    means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
    stds = synergia.bunch.Core_diagnostics().calculate_std(bunch, means)
    print >>fo, "%20s   %20s   %20s"%("coord","mean","rms")
    print >>fo, "%20s   %20s   %20s"%("====================",
                                      "====================",
                                      "====================")
    for i in range(6):
        print >>fo, "%20s   %20.12e   %20.12e"%(coord_names[i], means[i], stds[i]\
)


################################################################################

def get_lattice():
    # read the lattice in from a MadX sequence file
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "cfoborodobo32.madx")
    lattice.set_all_string_attribute("extractor_type", "libff")
    return lattice

################################################################################

# enable or disable RF.
# lattice_elements is an iteratable collection of lattice elements
def setup_rf(enable, elements, rf_volt, harmon, lag, logger):

    ncavities = 0
    for elem in elements:
        if elem.get_type() == ET.rfcavity:
            ncavities = ncavities + 1
            if enable:
                elem.set_double_attribute("volt", rf_volt)
                elem.set_double_attribute("lag", lag)
                elem.set_double_attribute("harmon", harmon)
                # frequency will be set by tune lattice
            else:
                elem.set_double_attribute("volt")
                elem.set_double_attribute("lag", 0.0)
                elem.set_double_attribute("harmon", 0.0)
    print("Set parameters for ", ncavities, " RF cavities")

    for elem in elements:
        if elem.get_type() == ET.rfcavity:
            print("rfcavity: ", elem, file=logger)

    # print('closed orbit: ', closed_orbit, file=logger)

    return

################################################################################

def main():

    logger = synergia.utils.Logger(0)
    lattice = get_lattice()
    print('Read lattice, length = {}, {} elements'.format(lattice.get_length(), len(lattice.get_elements())), file=logger)

    # set the momentum of the beam.  This could have been in a beam statement
    # in the lattice file, but this gives the option of setting it on the
    # command line by creating a reference particle with the desired momentum.

    # create the Reference particle object with the correct energy
    total_energy = opts.energy + synergia.foundation.pconstants.mp

    refpart = synergia.foundation.Reference_particle(1, synergia.foundation.pconstants.mp, total_energy)

    # set it into the lattice object
    lattice.set_reference_particle(refpart)

    energy = refpart.get_total_energy()
    momentum = refpart.get_momentum()
    gamma = refpart.get_gamma()
    beta = refpart.get_beta()

    print("energy: ", energy, file=logger)
    print("momentum: ", momentum, file=logger)
    print("gamma: ", gamma, file=logger)
    print("beta: ", beta, file=logger)

    # First set up RF without acceleration to get closed orbit and
    # initialize bunch distribution
    if opts.enable_rf:
        setup_rf(True, lattice.get_elements(), opts.rf_volt, opts.harmon, 0, logger)
        closed_orbit = synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
        orbit_length = closed_orbit[4] * beta
    else:
        orbit_length = lattice.get_length()

    print('Orbit length: ', orbit_length, file=logger)

    f = open("cfoborodobo32_lattice.out", "w")
    print(lattice, file=f)
    f.close()

    num_bunches = opts.num_bunches
    bucket_length = orbit_length/opts.harmon
    
    macroparticles = opts.macroparticles
    real_particles = opts.real_particles
    print("macroparticles: ", macroparticles, file=logger)
    print("real_particles: ", real_particles, file=logger)

    sim = synergia.simulation.Bunch_simulator.create_bunch_train_simulator(
        refpart, macroparticles, real_particles, num_bunches, bucket_length)

    if opts.periodic:
        sim.set_longitudinal_boundary(
                synergia.bunch.LongitudinalBoundary.periodic, bucket_length)

    comm = synergia.utils.parallel_utils.Commxx()

    map = SIM.Lattice_simulator.get_linear_one_turn_map(lattice)

    print('map:', file=logger)
    print(np.array2string(map, max_line_width=200), file=logger)

    [l, v] = np.linalg.eig(map)
    print("eigenvalues: ", file=logger)
    for z in l:
        print("|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi), file=logger)

    [ax, bx, qx] = map2twiss(map[0:2,0:2])
    [ay, by, qy] = map2twiss(map[2:4, 2:4])
    [az, bz, qz] = map2twiss(map[4:6,4:6])

    print("Lattice parameters (assuming uncoupled map)", file=logger)
    print("alpha_x: ", ax, " alpha_y: ", ay, file=logger)
    print("beta_x: ", bx, " beta_y: ", by, file=logger)
    print("q_x: ", qx, " q_y: ", qy, file=logger)
    print("q_z: ", qz, " beta_z: ", bz, file=logger)

    (orig_xtune, orig_ytune, orig_cdt) = SIM.Lattice_simulator.calculate_tune_and_cdt(lattice)
    print("Original base tunes, x: ", orig_xtune, " y: ", orig_ytune, file=logger)

    do_adjust_tunes = False
    if opts.xtune or opts.ytune:
        do_adjust_tunes = True
        nh, nv = mark_fd_quads(lattice)
        print(nh, ' horizontal correctors')
        print(nv, ' vertical correctors')

        if opts.xtune:
            target_xtune = opts.xtune
        else:
            target_xtune = orig_xtune
        if opts.ytune:
            target_ytune = opts.ytune
        else:
            target_ytune = orig_ytune

    if do_adjust_tunes:
        print("adjusting tunes, x: ", target_xtune," y: ", target_ytune, file=logger)
        SIM.Lattice_simulator.adjust_tunes(lattice, target_xtune, target_ytune, 1.0e-6)
        (new_xtune, new_ytune, new_cdt) = SIM.Lattice_simulator.calculate_tune_and_cdt(lattice)
        print("Adjusted tunes, x: ", new_xtune, " y: ", new_ytune, file=logger)
        

    chrom = SIM.Lattice_simulator.get_chromaticities(lattice)
    hchrom = chrom.horizontal_chromaticity
    vchrom = chrom.vertical_chromaticity
    print('initial horizontal chromaticity: ', hchrom, file=logger)
    print('initial vertical chromaticity: ', vchrom, file=logger)

    adjust_chromaticity = False
    if opts.set_hchrom or opts.set_vchrom:
        nh, nv = mark_fd_sext(lattice)
        print('number sextupole correctors: ', nh, nv, file=logger)
        if not nh and not nv:
            raise RuntimeError('No sextupole correctors available')
        print('initial correctors', file=logger)
        print_fd_sext(lattice, logger)
        adjust_chromaticity = True
        if opts.set_hchrom:
            target_hchrom = opts.set_hchrom
        else:
            target_hchrom = hchrom
        if opts.set_vchrom:
            target_vchrom = opts.set_vchrom
        else:
            target_vchrom = vchrom

        print('adjusting chromaticities to: h: ', target_hchrom, ', v: ', target_vchrom, file=logger)
        SIM.Lattice_simulator.adjust_chromaticities(lattice, target_hchrom, target_vchrom)
        # read back new chromaticity
        chrom = SIM.Lattice_simulator.get_chromaticities(lattice)
        hchrom = chrom.horizontal_chromaticity
        vchrom = chrom.vertical_chromaticity
        print('final correctors', file=logger)
        print_fd_sext(lattice, logger)

    alpha_c = chrom.momentum_compaction
    slip_factor = chrom.slip_factor
    print('final horizontal chromaticity: ', hchrom, file=logger)
    print('final vertical chromaticity: ', vchrom, file=logger)
    print("alpha_c: ", alpha_c, ", slip_factor: ", slip_factor, file=logger)

    dist = synergia.foundation.Random_distribution(opts.seed, comm) 

    stdx = opts.stdx
    stdy = opts.stdy
    stdz = opts.stdz
    stddpop = opts.stddpop

    print("stdx: ", stdx, file=logger)
    print("stdy: ", stdy, file=logger)
    print("stdz: ", stdz, file=logger)
    print("stddpop: ", stddpop, file=logger)

    # populate bunches a 6D matched bunch using either normal forms or a 6D moments procedure

    if opts.matching == "6dmoments":
        print("Matching with 6d moments", file=logger)

        covars = synergia.bunch.get_correlation_matrix(map, stdx, stdy, stdz, beta, (0,2,4))
        means = np.zeros(6, dtype='d')
        print(file=logger)
        print('covariance matrix:', file=logger)
        print(np.array2string(covars), file=logger, flush=True)
        for b in range(num_bunches):
            bunch = sim.get_bunch(0, b)
            synergia.bunch.populate_6d(dist, bunch, means, covars)
            print_statistics(bunch, logger)

    elif opts.matching == "uniform":
        print("Transversely matched, longitudinally uniform beam", file=logger)
        covars = synergia.bunch.get_correlation_matrix(map, stdx, stdy, stddpop, beta, (0,2,5))
        means = np.zeros(6, dtype='d')
        print(file=logger)
        print(np.array2string(covars), file=logger, flush=True)
        for b in range(num_bunches):
            bunch = sim.get_bunch(0, b)
            synergia.bunch.populate_transverse_gaussian(dist, bunch, means, covars, bucket_length/beta)
            print_statistics(bunch, logger)
        
    else:
        # no other matching options for now
        pass

    gridx = opts.gridx
    gridy = opts.gridy
    gridz = opts.gridz

    steps = opts.steps
    turns = opts.turns

    grid = [gridx, gridy, gridz]

    # space charge
    if opts.spacecharge:
        # space charge active
        comm_group_size = opts.comm_group_size

        if ops.spacecharge == "3dopen-hockney":
            sc_ops = synergia.collective.Space_charge_3d_open_hockney_options(gridx, gridy, gridz)
            sc_ops.comm_group_size = comm_group_size

        elif solver == "2dopen-hockney":
            sc_ops = synergia.collective.Space_charge_2d_open_hockney_options(gridx, gridy, gridz)
            sc_ops.comm_group_size = comm_group_size

        #elif solver == "2dbassetti-erskine":
        #    space_charge = synergia.collective.Space_charge_2d_bassetti_erskine()

        elif solver == "rectangular":
            space_charge = synergia.collective.Space_charge_rectangular(grid, opts.pipe_size)
            sc_ops.comm_group_size = comm_group_size

        else:
            sys.stderr.write("foborodobo32.py: solver must be either 3dopen-hockney, 2dopen-hockney, or rectangular\n")
            sys.exit(1)

        stepper = synergia.simulation.Split_operator_stepper(sc_ops, steps)

    else:
        # space charge not active

        if opts.stepper == "splitoperator":
            sc_ops = synergia.collective.Dummy_CO_options()
            stepper = synergia.simulation.Split_operator_stepper(sc_ops, steps)

        #elif opts.stepper == "independent":
        #    stepper = synergia.simulation.Independent_stepper(steps)

        elif opts.stepper == "elements":
            stepper = synergia.simulation.Independent_stepper_elements(steps)

        else:
            sys.stderr.write("mi.py: stepper must be either splitopertor,independent, or elements\n")
            sys.exit(1)

    # propagator
    propagator = synergia.simulation.Propagator(lattice, stepper)

    # diagnostics for the bunches
    for bunch_num in range(num_bunches):

        # diagnostics always on
        sim.reg_diag_per_turn(synergia.bunch.Diagnostics_full2("diag_b%03d.h5"%bunch_num),
                              bunch_idx = bunch_num)

        if opts.step_basic:
            sim.reg_diag_per_step(synergia.bunch.Diagnostics_basic("mi_step_basic_b%03d.h5"%bunch_num),
                    bunch_idx = bunch_num)

        if opts.particles:
            if opts.particles_period == 0:
                sim.reg_diag_per_turn(synergia.bunch.Diagnostics_particles("particles_b%03d.h5"%bunch_num), bunch_idx = bunch_num)
            else:
                turn_list = list(range(0, turns, opts.particles_period))
                if turns-1 not in turn_list:
                    turn_list.append(turns-1)
                    sim.reg_diag_turn_listed(synergia.bunch.Diagnostics_particles("mi_particles_b%03d.h5"%bunch_num), 
                                             bunch_idx = bunch_num, turns = turn_list)

        # enable track saving
        # each processor will save tracks/proc tracks
        if opts.tracks:
            trkfile = 'tracks_b%03d.h5'%bunch_num
            sim.reg_diag_per_turn(synergia.bunch.Diagnostics_bulk_track(trkfile, opts.tracks), bunch_idx = bunch_num)

    # max simulation turns
    sim.set_max_turns(opts.turns)

    # logs
    #simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_STEP)
    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_TURN)
    screen = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.DEBUG)

    print('number of lattice elements: ', len(lattice.get_elements()))
    print('lattice first element', lattice.get_elements()[0])
    print()
    lattice.get_elements()[0].set_double_attribute('a1', -9999.0)
    print('lattice first element mod: ', lattice.get_elements()[0])
    print()
    print()
    print()
    print('len(propagator.get_lattice_mutable().get_elements()): ', len(propagator.get_lattice_mutable().get_elements()))


    lattice_mutable = propagator.get_lattice_mutable()
    print('number propagator.get_lattice_mutable(): ', len(lattice_mutable.get_elements()))
    lattice_const = propagator.get_lattice_const()
    print('number of propagator.get_lattice_const(): ', len(lattice_const.get_elements()))
    print('number of propagator lattice elements: ', len(propagator.get_lattice_elements()))
    print()
    print()
    print('dir(lattice_mutable): ', dir(lattice_mutable))
    print()
    print('dir(lattice_const): ', dir(lattice_const))
    print()
    print('propagator lattice_elements first element orig: ', propagator.get_lattice_elements()[0])
    print()
    propagator.get_lattice_elements()[0].set_double_attribute('a1', 9999.0)
    print('propagator lattice_elements first element mod: ', propagator.get_lattice_elements()[0])
    print()
    print()
    print()
    print('lattice_mutable first element orig: ', lattice_mutable.get_elements()[0])
    lattice_mutable.get_elements()[0].set_double_attribute('a1', -222.0)
    print('lattice mutable first element mod: ', lattice_mutable.get_elements()[0])
    print()
    print('lattice_const first element orig: ', lattice_const.get_elements()[0])
    lattice_const.get_elements()[0].set_double_attribute('a1', 333.0)
    print('lattice const first element mod: ', lattice_const.get_elements()[0])
    print()
    print()
    print('propagator slices')
    i=0
    slices = propagator.get_lattice_element_slices()
    for s in slices:
        if i >= 10:
            break
        print(s, s)
        i = i+1

    # activate acceleration
    # setup_rf(True, propagator.get_lattice_elements(), opts.rf_volt, opts.harmon, opts.lag, logger)

    # propagate
    # propagator.propagate(sim, simlog, turns)

if __name__ == "__main__":
    main()
    #try:
    #    main()
    #except Exception as e:
    #    sys.stderr.write(str(e) + '\n')
    #    synergia.MPI.COMM_WORLD.Abort(777)


