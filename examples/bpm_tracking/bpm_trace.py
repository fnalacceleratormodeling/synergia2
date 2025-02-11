#!/usr/bin/env python

import sys, os
import numpy as np
import synergia
import synergia.simulation as SIM

import mpi4py.MPI as MPI

ET =  synergia.lattice.element_type

lattice_json = "booster-00000.json"
macroparticles = 1000
realparticles = 4.5e12/81
emitx = 7.0e-6 # normalized x emittance 95%
emity = 7.0e-6 # normalized y emittance 95%
stddE = 2.0e-3 # 4 sigma dE standard deviation
rndseed = 12345679
turns = 1000

def get_lattice():
    with open(lattice_json, 'rb') as f:
        lattice = synergia.lattice.Lattice.load_from_json(f.read())
    return lattice

def print_lattice(lattice):

    print('lattice length: ', lattice.get_length(), len(lattice.get_elements()), 'elements')

    refpart = lattice.get_reference_particle()
    print('reference particle characteristics')
    print('energy: ', refpart.get_total_energy())
    print('momentum: ', refpart.get_momentum())
    print('gamma: ', refpart.get_gamma())
    print('beta: ', refpart.get_beta())
    print('mass: ', refpart.get_mass())

    print()
    print('RF cavities:')
    s = 0.0

    for elem in lattice.get_elements():
        if elem.get_type() == ET.rfcavity:
            nm = elem.get_name()
            V = elem.get_double_attribute('volt')*1.0e-3  # GV
            freq = elem.get_double_attribute('freq')*1.0e6 # Hz
            phase = elem.get_double_attribute('lag') * 2.0*np.pi # radians
            print(f'    {s:12.9f} {nm} V:{V:8.3g}GV  freq:{freq}Hz phase:{phase}')
        s = s + elem.get_length()

    print()
    print('Sextupoles')

    s = 0.0
    for elem in lattice.get_elements():
        if elem.get_type() == ET.sextupole:
            nm = elem.get_name()
            k2 = elem.get_double_attribute('k2')
            print(f'    {s:12.9f} {nm} k2:{k2}')
        s = s + elem.get_length()

    # save lattice functions
    SIM.Lattice_simulator.CourantSnyderLatticeFunctions(lattice)
    SIM.Lattice_simulator.calc_dispersions(lattice)
            
    f = open('lattice_functions.csv', 'w')
    print('# name s beta_x alpha_x psi_x Dx D\'x beta_y alpha_y psi_y', file=f)
    
    for elem in lattice.get_elements():
        print(elem.lf.arcLength, " ", end='', file=f)
        print(elem.lf.beta.hor, " ", end='', file=f)
        print(elem.lf.alpha.hor, " ", end='', file=f)
        print(elem.lf.psi.hor, " ", end='', file=f)
        print(elem.lf.dispersion.hor, " ", end='', file=f)
        print(elem.lf.dPrime.hor, " ", end='', file=f)
        print(elem.lf.beta.ver, " ", end='', file=f)
        print(elem.lf.alpha.ver, " ", end='', file=f)
        print(elem.lf.psi.ver, file=f)
        pass

    f.close()
    return

def create_and_populate(lattice, emitx, emity, stddE):

    refpart = lattice.get_reference_particle()
    sim = SIM.Bunch_simulator.create_single_bunch_simulator(
        refpart, macroparticles, realparticles)

    bunch = sim.get_bunch(0, 0)

    momentum = refpart.get_momentum()
    beta =  refpart.get_beta()
    betagamma = beta*refpart.get_gamma()
    beta_x = lattice.get_elements()[-1].lf.beta.hor
    Dx = lattice.get_elements()[-1].lf.dispersion.hor
    beta_y = lattice.get_elements()[-1].lf.beta.ver
    
    stddpop = stddE/(4 * beta * momentum)
    stdx = np.sqrt((beta_x * emitx)/(4*betagamma) + Dx**2 * stddpop**2)
    stdy = np.sqrt(beta_y * emity/(4*betagamma))

    map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)
    covars = synergia.bunch.get_correlation_matrix(map, stdx, stdy, stddpop, beta, (0,2,5))
    means = np.zeros(6, dtype='d')
    dist = synergia.foundation.Random_distribution(rndseed, MPI.COMM_WORLD.rank)
    synergia.bunch.populate_6d(dist, bunch, means, covars)

    # start with a kick
    parts = bunch.get_particles_numpy()
    parts[:, 1] = 1.0e-4
    parts[:, 3] = -1.0e-4

    return sim

def print_bunch_statistics(bunch, logger=sys.stdout):

    parts = bunch.get_particles_numpy()
    print(parts.shape,  ", ", parts.size , file=logger)
    print("shape: {0}, {1}".format(parts.shape[0], parts.shape[1]), file=logger)

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch, mean)
    print("mean = {}".format(mean), file=logger)
    print("std = {}".format(std), file=logger)
    return

# register diagnostics to save mean position at each BPM.
# In this particular lattice all the BPMs have the same name. The
# reg_diag_at_element is triggered by the element name so all the BPM
# data will end up in the same file but on the good side I only have
# to add diagnostics one time for all the BPMs.
def register_diagnostics(sim, lattice):
    for elem in lattice.get_elements():
        nm = elem.get_name()
        if nm == "bpms":
            diag = synergia.bunch.Diagnostics_full2(f'bpmdata.h5')
            sim.reg_diag_at_element(diag, elem)
            break
                
    diag_full = synergia.bunch.Diagnostics_full2("diag.h5")
    sim.reg_diag_per_turn(diag_full)

    return


def get_propagator(lattice):
    stepper = SIM.Independent_stepper_elements(1)
    propagator = SIM.Propagator(lattice, stepper)
    return propagator

def main():
    lattice = get_lattice()
    refpart = lattice.get_reference_particle()

    print_lattice(lattice)

    sim = create_and_populate(lattice, emitx, emity, stddE)
    
    print_bunch_statistics(sim.get_bunch(0, 0), sys.stdout)

    register_diagnostics(sim, lattice)

    propagator = get_propagator(lattice)
    
    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_TURN)

    propagator.propagate(sim, simlog, turns)

    return

if __name__ == "__main__":
    main()
