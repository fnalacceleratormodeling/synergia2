
#import sys

from mpi4py import MPI
import numpy as np
import synergia

from mi_multibunch_options import opts

def print_statistics(bunch):

    parts = bunch.get_particles_numpy()
    print(parts.shape,  ", ", parts.size )
    print("shape: {0}, {1}".format(parts.shape[0], parts.shape[1]))

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch, mean)
    print("mean = {}".format(mean))
    print("std = {}".format(std))

def get_lattice_mi():
    lsexpr = synergia.utils.pylsexpr.read_lsexpr_file("mi20_ra_08182020_tuned.lsx")
    lattice = synergia.lattice.Lattice(lsexpr)
    lattice.set_all_string_attribute("extractor_type", "libff")
    #synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    return lattice

def run_mi():

    # parse options
    #opts.parse_argv(sys.argv)

    solver = opts.solver
    real_particles = opts.real_particles
    gridx = opts.gridx
    gridy = opts.gridy
    gridz = opts.gridz
    partpercell = opts.partpercell
    num_bunches = opts.num_bunches
    steps = opts.steps
    turns = opts.turns

    grid = [gridx, gridy, gridz]
    macro_particles = gridx * gridy * gridz * partpercell
    spacing = opts.bucket_length

    seed = 4

    # lattice
    lattice = get_lattice_mi()
    ref_part = lattice.get_reference_particle()
    #print(lattice)

    # bunch simulator
    sim = synergia.simulation.Bunch_simulator.create_bunch_train_simulator(
            ref_part, macro_particles, real_particles, num_bunches, spacing)

    # periodic
    if opts.periodic:
        sim.set_longitudinal_boundary(
                synergia.bunch.LongitudinalBoundary.periodic, spacing)

    # populate particles
    means = np.zeros([6], 'd')
    covars = np.load("correlation_matrix.npy")[()]
    #print(covars)

    comm = synergia.utils.parallel_utils.Commxx()
    dist = synergia.foundation.Random_distribution(seed, comm) 

    for b in range(num_bunches):
        bunch = sim.get_bunch(0, b)
        synergia.bunch.populate_6d(dist, bunch, means, covars)
        print_statistics(bunch)

    # space charge
    if opts.spacecharge:

        comm_group_size = opts.comm_group_size

        if solver == "3dopen-hockney":
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
            sys.stderr.write("mi.py: solver must be either 3dopen-hockney, 2dopen-hockney, or rectangular\n")
            sys.exit(1)

        stepper = synergia.simulation.Split_operator_stepper(sc_ops, steps)

    else:

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

    # diagnostics
    for bunch_num in range(num_bunches):
        #if opts.lam_dump:
            #bunch_train_simulator.add_per_forced_diagnostics_step(bunch_num, 
            #        synergia.bunch.Diagnostics_particles("lam_particles_b%03d.h5"%bunch_num), 
            #        opts.lam_period)

        for elem in lattice.get_elements():
            if elem.get_string_attribute("force_diagnostics", "") == "true":
                sim.reg_diag_per_element(
                        synergia.bunch.Diagnostics_particles("lam_particles_b%03d.h5"%bunch_num), 
                        elem, opts.lam_period, 0, bunch_num)

        if opts.turn_full2:
            sim.reg_diag_per_turn(synergia.bunch.Diagnostics_full2("mi_full2_b%03d.h5"%bunch_num),
                    0, bunch_num)

        if opts.step_basic:
            sim.reg_diag_per_step(synergia.bunch.Diagnostics_basic("mi_step_basic_b%03d.h5"%bunch_num),
                    0, bunch_num)

        if opts.step_full2:
            sim.reg_diag_per_step(synergia.bunch.Diagnostics_full2("mi_step_full2_b%03d.h5"%bunch_num),
                    0, bunch_num)

        if opts.turn_particles:
            turn_list = list(range(0, turns, opts.particles_period))
            if turns-1 not in turn_list:
                turn_list.append(turns-1)
                sim.reg_diag_turn_listed(synergia.bunch.Diagnostics_particles("mi_particles_b%03d.h5"%bunch_num), 
                        0, bunch_num, turn_list)

        # enable track saving
        # each processor will save tracks/proc tracks
        if opts.turn_tracks:
            # save 2000 turns worth of tracking in chunks of 200 turns
            turns_per_file = 200
            total_turns = 2000
            for trkturn in range(0,total_turns,turns_per_file):
                turn0 = trkturn
                turn1 = turn0+turns_per_file
                trkfile = "tracks_b%03d-%d-%d.h5"%(bunch_num,turn0,turn1)
                sim.add_per_turn(synergia.bunch.Diagnostics_bulk_track(trkfile, opts.turn_tracks), 
                        0, bunch_num, list(range(turn0,turn1)))

    # max simulation turns
    sim.set_max_turns(opts.max_turns)

    # logs
    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_STEP)
    screen = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.DEBUG)

    # propagate
    propagator.propagate(sim, simlog, turns)

    # print simple timer
    synergia.utils.parallel_utils.simple_timer_print(screen)

def main():

    print("run sis_18")
    print("my rank =", MPI.COMM_WORLD.Get_rank())

    #run2()
    #checkpoint_resume()

    run_mi()

main()

