#!/usr/bin/env python
import sys
import os
import synergia
import mpi4py.MPI as MPI
import numpy as np

import basic_toolkit
import beamline
import mxyzptlk
import physics_toolkit

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


try:

    myrank = 0

    #============================================
    # read the lattice and
    # turn on the RF by setting voltage, frequency and phase in RF cavities

    harmno = 588
    lattice = synergia.lattice.Mad8_reader().get_lattice("ring_p_q605", "mi20-egs-thinrf.lat")
    lattice_length = lattice.get_length()
    if myrank == 0:
        print("original lattice length: ", lattice_length)


    reference_particle = lattice.get_reference_particle()
    energy = reference_particle.get_total_energy()
    beta = reference_particle.get_beta()
    gamma = reference_particle.get_gamma()

    brho = lattice.get_reference_particle().get_momentum()*1.0e9/synergia.foundation.pconstants.c

    if myrank == 0:
        print("lattice length: ", lattice_length)
        print("energy: ", energy)
        print("beta: ", beta)
        print("gamma: ", gamma)
        print("brho: ", brho)

    # set rf cavity frequency (units are MHz)
    # harmno * beta * c/ring_length
    freq = 1.0e-6 * harmno * beta * synergia.foundation.pconstants.c/lattice_length
    if myrank == 0:
        print("RF frequency: ", freq)


    if myrank == 0:
        print("Begin setting RF voltage...")

    # rf cavity voltage, is 1.0 MV total distributed over 18 cavities.  MAD8
    # expects cavities voltages in  units of MV.
    rf_voltage = 1.0/18

    for elem in lattice.get_elements():
        if elem.get_type() == "rfcavity":
            elem.set_double_attribute("volt", rf_voltage)
            elem.set_double_attribute("freq", freq)

    if myrank == 0:
        print("Finish setting RF voltage...")

    #============================================

    if myrank == 0:
        print("Adding multipole attributes to lattice elements:")

    for elem in lattice.get_elements():
        # use specific average values for elements
        ename = elem.get_name()
        etype = elem.get_type()

        if ename[0:3] == "iqb" and etype == "quadrupole":
            strength = elem.get_double_attribute("k1")
            if strength > 0.0:
                # set skew quad to avg value
                elem.set_double_attribute("a1", -6.45295832275e-05)
                elem.set_double_attribute("b2", -0.00478869825563)
                elem.set_double_attribute("a2", -0.00160487324701)
            else:
                elem.set_double_attribute("a1", 2.8334536489e-06)
                elem.set_double_attribute("b2", 0.00260279001331)
                elem.set_double_attribute("a2", -0.00146028752357)
        elif ename[0:3] == "iqc" and etype == "quadrupole":
            strength = elem.get_double_attribute("k1")
            if strength > 0.0:
                # set skew quad to avg value
                elem.set_double_attribute("a1", 6.060810826e-05)
                elem.set_double_attribute("b2", 0.000360356508614)
                elem.set_double_attribute("a2",  0.00473082479225)
            else:
                elem.set_double_attribute("a1", 8.67982864331e-05)
                elem.set_double_attribute("b2", 0.00182760395913)
                elem.set_double_attribute("a2", -0.00190608360211)
        elif ename[0:3] == "iqd" and etype == "quadrupole":
            strength = elem.get_double_attribute("k1")
            if strength > 0.0:
                # set skew quad to avg value
                elem.set_double_attribute("a1", -2.21172896738e-05)
                elem.set_double_attribute("b2", 0.00319614326359)
                elem.set_double_attribute("a2", -0.000546928339363)
            else:
                elem.set_double_attribute("a1", 0.00049969635545)
                elem.set_double_attribute("b2", -0.0021999403501)
                elem.set_double_attribute("a2", -0.00102807817247)
        elif ename[0:3] == "ida" and etype == "sbend":
            # set normal and skew quad to avg value
            elem.set_double_attribute("b1", -0.00250687841324)
            elem.set_double_attribute("a1", -0.00274228645289)
            elem.set_double_attribute("b2", -0.0960362907447)
            elem.set_double_attribute("a2", -0.0247651386382)
        elif ename[0:3] == "idb" and etype == "sbend":
            elem.set_double_attribute("b1", 0.00116758066456)
            elem.set_double_attribute("a1", 0.00885250490156)
            elem.set_double_attribute("b2", 0.0101144512846)
            elem.set_double_attribute("a2", 0.00885250490156)
        elif ename[0:3] == "idc" and etype == "sbend":
            elem.set_double_attribute("b1", -0.00239160727335)
            elem.set_double_attribute("a1", -0.00374635022688)
            elem.set_double_attribute("b2", -0.0133260815125)
            elem.set_double_attribute("a2", -0.0143108656016)
        elif ename[0:3] == "idd" and etype == "sbend":
            elem.set_double_attribute("b1", 0.000792891581156)
            elem.set_double_attribute("a1", -0.00149762985985)
            elem.set_double_attribute("b2", -0.129264349405)
            elem.set_double_attribute("a2", 0.0165762123704)

    # now handle corrector settings and magnet movements
    # the numbers all come from mi_orbit_bumps.py
    first_h522 = True
    corrector_length = 0.3048 # from mi_orbit_bumps.py
    for elem in lattice.get_elements():
        ename = elem.get_name()
        # corrector settings were specified in KGauss
        if ename == "h520":
            kick = corrector_length * (0.231776/10.0)/brho
            elem.set_double_attribute("kick", kick)
        elif ename == "h522" and first_h522:
            kick = corrector_length * (-0.026296630050395536/10.0)/brho
            elem.set_double_attribute("kick", kick)
            first_h522 = False
        elif ename == "h524":
            kick = corrector_length * (0.26285942292792525/10.0)/brho
            elem.set_double_attribute("kick", kick)
        # offsets are originally specified in mm in a left-handed system
        elif ename == "iqd016":
            elem.set_double_attribute("hoffset", -(-0.531/1000.0))
        elif ename == "iqc022":
            elem.set_double_attribute("hoffset", -(-2.377/1000.0))
        elif ename == "iqe072":
            elem.set_double_attribute("hoffset", -(-0.710/1000.0))
        elif ename == "iqd024":
            elem.set_double_attribute("hoffset", -(-2.398/1000.0))
        elif ename == "iqd018":
            elem.set_double_attribute("hoffset", -(-0.329/1000.0))

    # elements not to be simplified
    dont_simplify = ["lam52", "h520", "h522", "h524", "iqd016", "iqc022",
                     "iqe072", "iqd024", "iqd018"]

    for elem in lattice.get_elements():
        if elem.get_name() in dont_simplify:
            elem.set_string_attribute("no_simplify", "true")
        elem.set_string_attribute("extractor_type", "chef_propagate")

    lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 2)
    print("lattice is a ring: ", lattice_simulator.is_ring())
    
    synergia.lattice.xml_save_lattice(lattice, "fnal_main_injector.xml")

    map = lattice_simulator.get_linear_one_turn_map()
    if myrank==0:
    	print("one turn map from synergia2.5 infrastructure")
    	print(np.array2string(map, max_line_width=200))


    [l, v] = np.linalg.eig(map)

    if myrank==0:
    	print("eigenvalues: ")
    	for z in l:
            print("|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi))

    [ax, bx, qx] = map2twiss(map[0:2,0:2])
    [ay, by, qy] = map2twiss(map[2:4, 2:4])
    [az, bz, qz] = map2twiss(map[4:6,4:6])
    stdz=1.19
    dpop = stdz/bz

    if myrank == 0:
    	print("Lattice parameters (assuming uncoupled map)")
    	print("alpha_x: ", ax, " alpha_y: ", ay)
    	print("beta_x: ", bx, " beta_y: ", by)
    	print("q_x: ", qx, " q_y: ", qy)
    	print("beta_z: ", bz)
    	print("delta p/p: ", dpop)

    proton = beamline.Proton(energy)
    beamline = lattice_simulator.get_chef_lattice().get_beamline()

    print("Initial proton: ", proton.State())
    beamline.propagate(proton)
    print("Final proton: ", proton.State())

except Exception as e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
