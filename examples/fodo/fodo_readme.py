#!/usr/bin/env python3

import numpy as np
import synergia
from synergia.utils import Commxx

macroparticles = 1048576
real_particles = 2.94e12
turns = 10
gridx = 32
gridy = 32
gridz = 128

def get_lattice():

    fodo_madx = """
      beam, particle=proton,pc=3.0;

      o: drift, l=8.0;
      f: quadrupole, l=2.0, k1=0.071428571428571425;
      d: quadrupole, l=2.0, k1=-0.071428571428571425;

      fodo: sequence, l=20.0, refer=entry;
      fodo_1: f, at=0.0;
      fodo_2: o, at=2.0;
      fodo_3: d, at=10.0;
      fodo_4: o, at=12.0;
      endsequence;
    """

    reader = synergia.lattice.MadX_reader()
    reader.parse(fodo_madx)
    lattice = reader.get_lattice('fodo')
    return lattice

def create_simulator(ref_part):

    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            ref_part, macroparticles, real_particles)

    bunch = sim.get_bunch()

    bunch_means = np.zeros(6, dtype='d')
    bunch_covariances = np.array(
        [[3.0509743977035345e-05, 2.2014134466660509e-06, 0, 0, 0, 0],
         [2.2014134466660509e-06, 1.9161816525115869e-07, 0, 0, 0, 0],
         [0, 0, 7.5506914064526925e-06, -6.6846812465678249e-07, 0, 0],
         [0, 0, -6.6846812465678249e-07, 1.9161816525115867e-07, 0, 0],
         [0, 0, 0, 0, 0.00016427607645871527, 0],
         [0, 0, 0, 0, 0, 1e-08]])
    
    dist = synergia.foundation.PCG_random_distribution(1234567, Commxx.World.rank())

    synergia.bunch.populate_6d(dist, 
        bunch, 
        bunch_means,
        bunch_covariances)

    return sim

def create_propagator(lattice):
    sc_ops = synergia.collective.Space_charge_3d_open_hockney_options(gridx, gridy, gridz)
    sc_ops.comm_group_size = 1

    stepper = synergia.simulation.Split_operator_stepper_elements(sc_ops, 1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

def run():

    # logger
    screen = synergia.utils.parallel_utils.Logger(0, 
      synergia.utils.parallel_utils.LoggerV.INFO)

    simlog = synergia.utils.parallel_utils.Logger(0, 
      synergia.utils.parallel_utils.LoggerV.INFO_TURN)

    # components
    lattice = get_lattice()
    sim = create_simulator(lattice.get_reference_particle())
    propagator = create_propagator(lattice)

    # diagnostics
    diag_full2 = synergia.bunch.Diagnostics_full2("diag_full.h5")
    sim.reg_diag_per_turn(diag_full2)

    # propagate
    propagator.propagate(sim, simlog, turns)

    # save
    synergia.simulation.checkpoint_save(propagator, sim)

def main():

    try:
        run()
    except:
        raise RuntimeError("Failure to launch fodo.run")

main()
