from mpi4py import MPI
import numpy as np
import synergia

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
    print p.shape

    #diag_full2 = synergia.bunch.Diagnostics_full2("diag_full.h5")
    #sim.reg_diag_per_turn("full2", diag_full2)

    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO)
    propagator.propagate(sim, simlog)


def main():

    print("run sis_18")
    print "my rank =", MPI.COMM_WORLD.Get_rank()
    run()

main()

