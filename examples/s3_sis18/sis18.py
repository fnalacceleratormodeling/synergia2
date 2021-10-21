from mpi4py import MPI
import numpy as np
import synergia

class mydiag(synergia.bunch.Diagnostics):
    def __init__(self, filename):
        synergia.bunch.Diagnostics.__init__(self, "mydiag", filename)
        print("mydiag created")

    def do_update(self, bunch):
        print("my diag update")

    def do_reduce(self, comm, root):
        print("my diag reduce")


def print_statistics(bunch):

    parts = bunch.get_particles_numpy()
    print(parts.shape,  ", ", parts.size )
    print("shape: {0}, {1}".format(parts.shape[0], parts.shape[1]))

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch, mean)
    print("mean = {}".format(mean))
    print("std = {}".format(std))


def get_lattice():

    #lsexpr = synergia.utils.pylsexpr.read_lsexpr_file("sis18-6.lsx")
    #lattice = synergia.lattice.Lattice(lsexpr)

    reader = synergia.lattice.MadX_reader()
    reader.parse_file("sis18.madx")
    lattice = reader.get_lattice("machine")

    lattice.set_all_string_attribute("extractor_type", "libff")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    lattice.export_madx_file("exported_lattice.madx")

    return lattice

def create_simulator(ref_part):
    #comm = synergia.utils.parallel_utils.Commxx()
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            ref_part, 4194394, 2.94e10)

    bunch = sim.get_bunch()
    #bunch.read_file("bunch_particles_4M.h5")

    means = np.zeros([6], 'd') 
    covars = np.load("correlation_matrix.npy")[()]

    print("means: ", means)
    print("covars: ", covars)

    comm = synergia.utils.parallel_utils.Commxx()
    #dist = synergia.foundation.Random_distribution(5, comm)
    dist = synergia.foundation.PCG_random_distribution(5, comm)

    synergia.bunch.populate_6d(dist, bunch, means, covars)

    return sim

def create_propagator(lattice):
    sc_ops = synergia.collective.Space_charge_2d_open_hockney_options(64, 64, 64)
    sc_ops.comm_group_size = 1

    stepper = synergia.simulation.Split_operator_stepper(sc_ops, 71)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

def run2():

    screen = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.DEBUG)

    lattice = get_lattice()
    sim = create_simulator(lattice.get_reference_particle())
    propagator = create_propagator(lattice)

    sim.get_bunch().print_statistics(screen);

    class context:
        steps = 0

    def action(sim, lattice, turn, step):
        #nonlocal another_steps
        #another_steps += 1
        context.steps += 1

    sim.reg_prop_action_step_end(action)

    # diagnostics
    diag_full2 = synergia.bunch.Diagnostics_full2("diag_full_py.h5")
    sim.reg_diag_per_turn(diag_full2)

    diag_bt = synergia.bunch.Diagnostics_bulk_track("diag_bt_py.h5", 1000, 0)
    sim.reg_diag_per_turn(diag_bt)

    diag_part = synergia.bunch.Diagnostics_particles("diag_part_py.h5", 100)
    sim.reg_diag_per_turn(diag_part)

    sim.reg_diag_per_turn(mydiag("mydiag.h5"))
    #diag_dummy = synergia.bunch.Diagnostics_dummy()
    #sim.reg_diag_per_turn(diag_dummy)

    # logger
    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_STEP)

    # propagate
    propagator.propagate(sim, simlog, 1)

    print("total steps = ", context.steps)
    sim.get_bunch().print_statistics(screen);

    # save
    synergia.simulation.checkpoint_save(propagator, sim);

def checkpoint_resume():

    [propagator, sim] = synergia.simulation.checkpoint_load();

    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_STEP)

    screen = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.DEBUG)

    sim.get_bunch().print_statistics(screen);

    propagator.propagate(sim, simlog, 1)

    sim.get_bunch().print_statistics(screen);
    synergia.utils.parallel_utils.simple_timer_print(screen)

def run():
    lsexpr = synergia.utils.pylsexpr.read_lsexpr_file("sis18-6.lsx")
    lattice = synergia.lattice.Lattice(lsexpr)
    lattice.set_all_string_attribute("extractor_type", "libff")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    ref = lattice.get_reference_particle()

    #print(lattice)

    sc_ops = synergia.collective.Space_charge_2d_open_hockney_options(64, 64, 64)
    sc_ops.comm_group_size = 1

    stepper = synergia.simulation.Split_operator_stepper(sc_ops, 71)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    #comm = synergia.utils.parallel_utils.Commxx()
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            ref, 4194394, 2.94e10)
    sim.set_turns(0, 1)

    bunch = sim.get_bunch()
    bunch.read_file("turn_particles_0000_4M.h5")


    parts = bunch.get_host_particles()
    p = np.array(parts, copy=False)
    print(p.shape,  ", ", p.size )
    print("shape: {0}, {1}".format(p.shape[0], p.shape[1]))

    print_statistics(bunch)

    class context:
        steps = 0

    def action(sim, lattice, turn, step):
        #nonlocal steps
        #steps += 1
        context.steps += 1
        print(context.steps)


    sim.reg_prop_action_step_end(action)

    diag_full2 = synergia.bunch.Diagnostics_full2("diag_full.h5")
    sim.reg_diag_per_turn(diag_full2)

    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO)
    propagator.propagate(sim, simlog)

    print(context.steps)
    print_statistics(bunch)

def main():

    print("run sis_18")
    print("my rank =", MPI.COMM_WORLD.Get_rank())
    run2()
    checkpoint_resume()

main()

