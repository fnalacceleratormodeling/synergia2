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


def get_lattice():
    lsexpr = synergia.utils.pylsexpr.read_lsexpr_file("sis18-6.lsx")
    lattice = synergia.lattice.Lattice(lsexpr)
    lattice.set_all_string_attribute("extractor_type", "libff")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    return lattice

def create_simulator(ref_part):
    #comm = synergia.utils.parallel_utils.Commxx()
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            ref_part, 4194394, 2.94e10)

    bunch = sim.get_bunch()
    bunch.read_file("turn_particles_0000_4M.h5")

    return sim

def create_propagator(lattice):
    sc_ops = synergia.collective.Space_charge_2d_open_hockney_options(64, 64, 64)
    sc_ops.comm_group_size = 1

    stepper = synergia.simulation.Split_operator_stepper(sc_ops, 71)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

def run2():
    lattice = get_lattice()
    sim = create_simulator(lattice.get_reference_particle())
    propagator = create_propagator(lattice)

    print_statistics(sim.get_bunch())

    class context:
        steps = 0

    def action(sim, lattice, turn, step):
        #nonlocal another_steps
        #another_steps += 1
        context.steps += 1

    sim.reg_prop_action_step_end(action)
    sim.set_turns(0, 1)

    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO)

    propagator.propagate(sim, simlog)

    print("total steps = ", context.steps)
    print_statistics(sim.get_bunch())

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

    #diag_full2 = synergia.bunch.Diagnostics_full2("diag_full.h5")
    #sim.reg_diag_per_turn("full2", diag_full2)

    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO)
    propagator.propagate(sim, simlog)

    print(context.steps)
    print_statistics(bunch)

def main():

    print("run sis_18")
    print("my rank =", MPI.COMM_WORLD.Get_rank())
    run2()

main()

