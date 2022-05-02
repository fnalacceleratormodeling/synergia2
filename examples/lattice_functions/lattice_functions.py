#!/usr/bin/env python3

import sys, os
import numpy as np
import synergia
import synergia.simulation as SIM
PCONST = synergia.foundation.pconstants

import matplotlib.pyplot as plt


def main():

    # read the lattice in from a MadX sequence file
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "foborodobo32.madx")
    print('Read lattice: ', lattice.get_name(), len(lattice.get_elements()), ' elements, length: ', lattice.get_length())
    #print(lattice)
    refpart = synergia.foundation.Reference_particle(1, PCONST.mp, 0.8 + PCONST.mp)
    lattice.set_reference_particle(refpart)
    print('energy: ', refpart.get_total_energy())
    print('momentum: ', refpart.get_momentum())
    print('gamma: ', refpart.get_gamma())
    print('beta: ', refpart.get_beta())
    print('mass: ', refpart.get_mass())

    SIM.Lattice_simulator.tune_circular_lattice(lattice)
    SIM.Lattice_simulator.CourantSnyderLatticeFunctions(lattice)
    SIM.Lattice_simulator.calc_dispersions(lattice)
    elem = lattice.get_elements()[0]
    # print('elem: ', dir(elem))
    # print('lf: ', dir(elem.lf))
    # print('lf.beta: ', dir(elem.lf.beta))
    # print('lf.arcLength: ', elem.lf.arcLength)
    # print('lf.dispersion: ', elem.lf.dispersion)
    arcLength = []
    beta_x = []
    alpha_x = []
    beta_y = []
    alpha_y = []
    psi_x = []
    psi_y = []
    disp_x = []
    disp_y = []
    for elem in lattice.get_elements():
        arcLength.append(elem.lf.arcLength)
        beta_x.append(elem.lf.beta.hor)
        alpha_x.append(elem.lf.alpha.hor)
        psi_x.append(elem.lf.psi.hor)
        disp_x.append(elem.lf.dispersion.hor)
        beta_y.append(elem.lf.beta.ver)
        alpha_y.append(elem.lf.alpha.ver)
        psi_y.append(elem.lf.psi.ver)
        disp_y.append(elem.lf.dispersion.ver)
    
    plt.figure()
    plt.title('beta x,y')
    plt.plot(arcLength, beta_x, label='beta_x')
    plt.plot(arcLength, beta_y, label='beta_y')
    plt.xlabel('s [m]')
    plt.ylabel('beta')
    plt.legend(loc='best')
    
    plt.figure()
    plt.title('alpha x,y')
    plt.plot(arcLength, alpha_x, label='alpha_x')
    plt.plot(arcLength, alpha_y, label='alpha_y')
    plt.xlabel('s [m]')
    plt.ylabel('alpha')
    plt.legend(loc='best')
    
    plt.figure()
    plt.title('phase x,y')
    plt.plot(arcLength, psi_x, label='psi_x')
    plt.plot(arcLength, psi_y, label='psi_y')
    plt.xlabel('s [m]')
    plt.ylabel('psi')
    plt.legend(loc='best')
    
    plt.figure()
    plt.title('dispersion x,y')
    plt.plot(arcLength, disp_x, label='dispersion x')
    plt.plot(arcLength, disp_y, label='dispersion y')
    plt.xlabel('s [m]')
    plt.ylabel('dispersion')
    plt.legend(loc='best')
    
    plt.show()
        
  
if __name__ == "__main__":
    main()
