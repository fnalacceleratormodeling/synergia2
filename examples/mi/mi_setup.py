#!/usr/bin/env python

import sys
import os
import synergia
import mpi4py.MPI as MPI
import numpy as np

import synergia.simulation as SIM 

#import mi_fixup

from mi_multibunch_options import opts

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

# get_fd_quads(lattice) reads input lattice and returns
# ( [list of focussing quad elements], [list of defocussing quad elements] )
def get_fd_quads(lattice):
    f_quads = []
    d_quads = []
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.quadrupole:
            ename = elem.get_name()
            if ((ename.find("iqb") == 0) or
                (ename.find("iqc") == 0) or
                (ename.find("iqd") == 0) or
                (ename.find("iqe") == 0) or
                (ename.find("iqf") == 0) or
                (ename.find("iqg") == 0)):
                if (elem.get_double_attribute("k1") > 0.0):
                    f_quads.append(elem)
                    elem.set_marker(synergia.lattice.marker_type.h_tunes_corrector)
                elif (elem.get_double_attribute("k1") < 0.0):
                    d_quads.append(elem)
                    elem.set_marker(synergia.lattice.marker_type.v_tunes_corrector)
    return (f_quads, d_quads)

######################################################

# get_fd_sextupoles(lattice) reads the input lattice nd returns
# ( [list of focussing sextupoles], [list of defocussing sextupoles] )
# for use in chromacity adjustment

def get_fd_sextupoles(lattice):
    f_sext = []
    d_sext = []
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.sextupole:
            if (elem.get_double_attribute("k2") > 0.0):
                f_sext.append(elem)
                elem.set_marker(synergia.lattice.marker_type.h_chrom_corrector)
            elif (elem.get_double_attribute("k2") < 0.0):
                d_sext.append(elem)
                elem.set_marker(synergia.lattice.marker_type.v_chrom_corrector)
    return (f_sext, d_sext)

#######################################################

def setup():

    try:
        logger = synergia.utils.parallel_utils.Logger(0, 
                synergia.utils.parallel_utils.LoggerV.DEBUG)

        #logger = synergia.utils.Logger(0)

        tjob_0 = MPI.Wtime()

        #memlogger = synergia.utils.Logger("memlog")


        #============================================
        # get parameters for run
        
        rf_voltage = opts.rf_voltage

        print("==== Run Summary ====", file=logger)
        print("RF Voltage: ", rf_voltage, file=logger)

        if opts.xtune_adjust:
            print("x tune adjustment: ", opts.xtune_adjust, file=logger)
        else:
            print("x tune adjustment: NONE", file=logger)
        if opts.ytune_adjust:
            print("y tune adjustment: ", opts.ytune_adjust, file=logger)
        else:
            print("y tune adjustment: NONE", file=logger)

        if opts.xchrom_adjust:
            print("x chromaticity adjustment: ", opts.xchrom_adjust, file=logger)
        else:
            print("x chromaticity adjustment: NONE", file=logger)

        if opts.ychrom_adjust:
            print("y chromaticity adjustment: ", opts.ychrom_adjust, file=logger)
        else:
            print("y chromaticity adjustment: NONE", file=logger)

        #============================================

        # read the lattice and
        # turn on the RF by setting voltage, frequency and phase in RF cavities

        """
        f = open('mi20_ra_08182020.lat')
        #f = open('testlattice.lat')
        lattice_lines = f.readlines()
        print("read ", len(lattice_lines), ' lines')
        mi_fixup.remove_comments_and_continuation(lattice_lines)

        mi_fixup.mi_fixup(lattice_lines)

        m8_reader = synergia.lattice.Mad8_reader()
        m8_reader.parse_string('\n'.join(lattice_lines))

        lattice = m8_reader.get_lattice('ring_p_q100')
        """

        lsexpr = synergia.utils.pylsexpr.read_lsexpr_file("mi20_raw.lsx")
        lattice = synergia.lattice.Lattice(lsexpr)
        lattice.set_all_string_attribute("extractor_type", "libff")

        print("lattice # elements: ", len(lattice.get_elements()))
        print("lattice length: ", lattice.get_length())

        harmno = 588

        if opts.start_element:
            print("reordering lattice to start at ", opts.start_element, file=logger)
            new_lattice = synergia.lattice.Lattice("mi", lattice.get_element_adaptor_sptr())
            copying = False
            start_i=-1
            elements = lattice.get_elements()
            for i in range(len(elements)):
                if elements[i].get_name() == opts.start_element:
                    start_i = i
                    break
            reorder_elements = elements[start_i:] + elements[:start_i]
            for elem in reorder_elements:
                new_lattice.append(elem)

            new_lattice.set_reference_particle(lattice.get_reference_particle())
            print("new lattice length: ", new_lattice.get_length(), file=logger)
            lattice = new_lattice

        #============================================

        # turn on the RF cavities
        reference_particle = lattice.get_reference_particle()
        energy = reference_particle.get_total_energy()
        beta = reference_particle.get_beta()
        gamma = reference_particle.get_gamma()

        print("energy: ", energy, file=logger)
        print("beta: ", beta, file=logger)
        print("gamma: ", gamma, file=logger)

        # set rf cavity frequency
        # harmno * beta * c/ring_length
        freq = harmno * beta * synergia.foundation.pconstants.c/lattice.get_length()
        # only for informational purposes
        print("RF frequency: ", freq, file=logger)

        print("Begin setting RF voltage...", file=logger)

        # rf cavity voltage, is 1.0 MV total distributed over 18 cavities.  MAD8
        # expects cavities voltages in  units of MV.
        for elem in lattice.get_elements():
            if elem.get_type() == synergia.lattice.element_type.rfcavity:
                elem.set_double_attribute("volt", rf_voltage)
                # set the harmonic number so the frequency is set 
                elem.set_double_attribute("harmon", harmno)
                # set the first pass frequency so I can get the bucket length
                elem.set_double_attribute("freq", freq*1.0e-6)

        print("Finish setting RF voltage...", file=logger)

        #============================================
        # adjust the tune and chromaticity of requested

        SIM.Lattice_simulator.set_closed_orbit_tolerance(1.0e-6)
        SIM.Lattice_simulator.tune_circular_lattice(lattice)

        tunes = SIM.Lattice_simulator.calculate_tune_and_cdt(lattice)
        xtune = tunes[0]
        ytune = tunes[1]

        chroms = SIM.Lattice_simulator.get_chromaticities(lattice)
        xchrom = chroms.horizontal_chromaticity
        ychrom = chroms.vertical_chromaticity


        # adjust the tunes
        print("Unadjusted x tune: ", xtune, file=logger)
        print("Unadjusted y tune: ", ytune, file=logger)

        if opts.xtune_adjust or opts.ytune_adjust:
            if opts.xtune_adjust:
                xtune += opts.xtune_adjust
            if opts.ytune_adjust:
                ytune += opts.ytune_adjust

            print("Adjusting tunes to:", file=logger)
            print("xtune: ", xtune, file=logger)
            print("ytune: ", ytune, file=logger)

            tune_tolerance = 1.0e-8

            f_quads, d_quads = get_fd_quads(lattice)
            print("There are ", len(f_quads), " focussing quadrupoles", file=logger)
            print("There are ", len(d_quads), " defocussing quadrupoles", file=logger)

            SIM.Lattice_simulator.adjust_tunes(
                    lattice, xtune, ytune, tune_tolerance)

            tunes = SIM.Lattice_simulator.calculate_tune_and_cdt(lattice)
            xtune = tunes[0]
            ytune = tunes[1]

            print("adjusted xtune: ", xtune, file=logger)
            print("adjusted ytune: ", ytune, file=logger)

        print("Unadjusted x chromaticity: ", xchrom, file=logger)
        print("Unadjusted y chromaticity: ", ychrom, file=logger)
            
        if opts.xchrom_adjust or opts.ychrom_adjust:

            if opts.xchrom_adjust:
                xchrom += opts.xchrom_adjust

            if opts.ychrom_adjust:
                ychrom += opts.ychrom_adjust

            print("Adjusting chromaticities to:", file=logger)
            print("x chrom: ", xchrom, file=logger)
            print("y chrom: ", ychrom, file=logger)

            f_sext, d_sext = get_fd_sextupoles(lattice)
            print("There are ", len(f_sext), " focussing sextupoles", file=logger)
            print("There are ", len(d_sext), " defocussing sextupoles", file=logger)
            
            SIM.Lattice_simulator.adjust_chromaticities(
                    lattice, xchrom, ychrom, 1.0e-6, 20)

            chroms = SIM.Lattice_simulator.get_chromaticities(lattice)
            xchrom = chroms.horizontal_chromaticity
            ychrom = chroms.vertical_chromaticity

            print("adjusted x chrom: ", xchrom, file=logger)
            print("adjusted y chrom: ", ychrom, file=logger)


        # The lattice is tuned now, write it out
        #synergia.utils.write_lsexpr_file(lattice.as_lsexpr(), "mi20_ra_08182020_tuned.lsx")

        # Get the covariance matrix
        map = SIM.Lattice_simulator.get_linear_one_turn_map(lattice)
        print("one turn map from synergia2.5 infrastructure", file=logger)
        print(np.array2string(map, max_line_width=200), file=logger)

        [l, v] = np.linalg.eig(map)

        #print "l: ", l
        #print "v: ", v

        print("eigenvalues: ", file=logger)
        for z in l:
            print("|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi), file=logger)

        [ax, bx, qx] = map2twiss(map[0:2,0:2])
        [ay, by, qy] = map2twiss(map[2:4, 2:4])
        [az, bz, qz] = map2twiss(map[4:6,4:6])
        stdz = opts.stdz
        dpop = stdz/bz

        print("Lattice parameters (assuming uncoupled map)", file=logger)
        print("alpha_x: ", ax, " alpha_y: ", ay, file=logger)
        print("beta_x: ", bx, " beta_y: ", by, file=logger)
        print("q_x: ", qx, " q_y: ", qy, file=logger)
        print("beta_z: ", bz, file=logger)
        print("delta p/p: ", dpop, file=logger)


        emitx = opts.norm_emit/(beta*gamma)
        emity = opts.norm_emit/(beta*gamma)

        stdx = np.sqrt(emitx*bx)
        stdy = np.sqrt(emity*by)

        print("emitx: ", emitx, file=logger)
        print("emity: ", emity, file=logger)
        print("target stdx: ", stdx, file=logger)
        print("target stdy: ", stdy, file=logger)

        #=========================================
        #
        # get the correlation matrix for bunch generation

        correlation_matrix = synergia.bunch.get_correlation_matrix(
                map, stdx, stdy, stdz/beta, beta, (0, 2, 4))

        np.save("correlation_matrix.npy", correlation_matrix)
        print(np.array2string(correlation_matrix, max_line_width=200), file=logger)

    except Exception as e:
        sys.stdout.write(str(e) + '\n')
        MPI.COMM_WORLD.Abort(777)

def main():

    print("running mi_setup")
    setup()

main()
