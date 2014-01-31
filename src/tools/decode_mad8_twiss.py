#!/usr/bin/env python

import sys
import os
import numpy as np
from fortranformat import FortranRecordReader as frr
import getopt
import synergia

DEBUG = False

def test_close(x, y, tolerance):
    if abs(x) < tolerance:
        return abs(x-y)<tolerance
    else:
        return abs((x-y)/x) < tolerance

def usage():
    print "usage: decode_madx_twiss [options] TWISS-file"
    print "    options: --xmlfile=file write xml description of lattice"
    print "             --lfcsvfile=file write lattice functions in csv format to file"
    print "             --lfnpfile=file write numpy save file of lattice functions"
    #print "             --lfcompare compare MAD8 lattice functions with Synergia lattice functions produced by the lattice"
    print "             --help  this message"

lattice_element_name_map = {
    'DRIF': 'drift', 'RBEN': 'rbend', 'SBEN': 'sbend', 'QUAD': 'quadrupole',
    'SEXT': 'sextupole', 'OCTU': 'octupole', 'MULT': 'multipole',
    'SOLE': 'solenoid', 'RFCA': 'rfcavity', 'ELSE': 'elseparator',
    'KICK': 'kicker', 'HKIC': 'hkicker', 'VKIC': 'vkicker', 'SROT': 'srot',
    'YROT': 'yrot', 'MONI': 'monitor', 'HMON': 'hmonitor', 'VMON': 'vmonitor',
    'MARK': 'marker', 'INST': 'instrument'}

def make_lattice_element(letype, lename, lelength, p1, p2, p3, p4, p5, p6, p7, p8):
    new_elem = synergia.lattice.Lattice_element(lattice_element_name_map[letype],
                                                lename.strip())
    new_elem.set_double_attribute('l', lelength)
    if (letype == 'DRIF') or (letype == "INST"):
        pass # nothing to do for drifts
    elif (letype == 'RBEN') or (letype == 'SBEN'):
        new_elem.set_double_attribute('angle', p1)
        if p2 != 0.0:
            new_elem.set_double_attribute('k1', p2)
        if p3 != 0.0:
            new_elem.set_double_attribute('k2', p3)
        if p4 != 0.0:
            new_elem.set_double_attribute('tilt', p4)
        if p5 != 0.0:
            new_elem.set_double_attribute('e1', p5)
        if p6 != 0.0:
            new_elem.set_double_attribute('e2', p6)
        if p7 != 0.0:
            new_elem.set_double_attribute('h1', p7)
        if p8 != 0.0:
            new_elem.set_double_attribute('h2', p8)
    elif letype == 'QUAD':
        if p2 != 0.0:
            new_elem.set_double_attribute('k1', p2)
        if p4 != 0.0:
            new_elem.set_double_attribute('tilt', p4)
    elif letype == 'SEXT':
        if p3 != 0.0:
            new_elem.set_double_attribute('k2', p3)
        if p4 != 0.0:
            new_elem.set_double_attribute('tilt', p4)
    elif letype == 'OCTU':
        if p4 != 0.0:
            new_elem.set_double_attribute('tilt', p4)
        if p5 != 0.0:
            new_elem.set_double_attribute('k3', p5)
    elif letype == 'MULT':
        new_elem.set_double_attribute('l', 0.0)
        if p1 != 0.0:
            new_elem.set_double_attribute('k0l', p1)
        if p2 != 0.0:
            new_elem.set_double_attribute('k1l', p2)
        if p3 != 0.0:
            new_elem.set_double_attribute('k2l', p3)
        if p4 != 0.0:
            new_elem.set_double_attribute('t0', p4)
        if p5 != 0.0:
            new_elem.set_double_attribute('k3l', p5)
        if p6 != 0.0:
            new_elem.set_double_attribute('t1', p6)
        if p7 != 0.0:
            new_elem.set_double_attribute('t2', p7)
        if p8 != 0.0:
            new_elem.set_double_attribute('t3', p8)
    elif letype == 'SOLE':
        new_elem.set_double_attribute('ks', p5)
    elif letype == 'RFCA':
        new_elem.set_double_attribute('freq', p5)
        new_elem.set_double_attribute('volt', p6)
        new_elem.set_double_attribute('lag', p7)
    elif letype == 'ELSE':
        new_elem.set_double_attribute('tilt', p4)
        new_elem.set_double_attribute('efield', p5)
    elif letype == 'KICK':
        new_elem.set_double_attribute('tilt', p4)
        new_elem.set_double_attribute('hkick', p5)
    elif (letype == 'HKIC') or (letype == 'VKIC'):
        new_elem.set_double_attribute('tilt', p4)
        new_elem.set_double_attribute('kick', p5)
    elif (letype == "SROT") or (letype == "YROT"):
        new_elem.set_double_attribute('l', 0.0)
        new_elem.set_double_attribute('angle', p5)
    elif (letype == "MONI") or (letype == "HMON") or (letype == "VMON"):
        pass # no additional attributes for monitors
    elif letype == "MARK":
        new_elem.set_double_attribute('l', 0.0)

    if DEBUG:
        print "new_elem: ", new_elem.as_string()
    return new_elem
    

# read twiss file returns a tuple with a lattice function and beamline information
def read_twiss_file(tfo):
    # read header
    hdrformat = frr('(5A8,I8,L8,A8)')
    progvrsn,datavrsn,madrundate, madruntime, madjobname,superperiod,symm,npos = hdrformat.read(tfo.readline())
    if DEBUG:
        print "progvrsn: ", progvrsn
        print "datavrsn: ", "->"+datavrsn+"<-"
        print "madrundate: ", "->"+madrundate+"<-"
        print "madruntime: ", "->"+madruntime+"<-"
        print "superperiod: ", superperiod
        print "symm: ", symm
        print "npos: ", "->"+npos+"<-"
    if datavrsn != "   TWISS":
        raise RuntimeError, "File is not a TWISS file"

    titleformat = frr('(A80)')
    runtitle = titleformat.read(tfo.readline())[0]
    if DEBUG:
        print "Run Title: ", runtitle
    num_records = int(npos)
    lfinfo = []
    lattice = synergia.lattice.Lattice()
    s = 0.0
    line1format = frr('(A4,A16,F12.6,3E16.9)')
    # keyword, name, , length, element specific
    line2format = frr('(5E16.9)')
    line3_5format = frr('(5E16.9)')
    # alf, bet, mu, dx, dpx, etc..
    for elem_num in range(num_records):
        letype, lename, lelength, p1, p2, p3 = line1format.read(tfo.readline())
        p4, p5, p6, p7, p8 = line2format.read(tfo.readline())
        if DEBUG:
            print letype, lename, lelength, p1, p2, p3
        
        alpha_x, beta_x, psi_x, D_x, Dprime_x = line3_5format.read(tfo.readline())
        alpha_y, beta_y, psi_y, D_y, Dprime_y = line3_5format.read(tfo.readline())
        x, px, y, py, suml = line3_5format.read(tfo.readline())
        # I use a dict because the simulation.Lattice_functions struct is
        # const
        lf = {}
        lf['name'] = lename
        lf['alpha_x'] = alpha_x
        lf['beta_x'] = beta_x
        lf['psi_x'] = psi_x
        lf['D_x'] = D_x
        lf['Dprime_x'] = Dprime_x
        lf['alpha_y'] = alpha_y
        lf['beta_y'] = beta_y
        lf['psi_y'] = psi_y
        lf['D_y'] = D_y
        lf['Dprime_y'] = Dprime_y
        lf['s'] = suml

        if elem_num != 0:
            lattice.append(make_lattice_element(letype, lename, lelength,
                                                p1, p2, p3, p4, p5, p6, p7, p8))
        lfinfo.append(lf)
        
    return lfinfo, lattice
    
def write_lfcsvfile(filename, lfinfo):
    lffo = open(filename, "w")
    print >>lffo, "#name s alpha_x beta_x psi_x alpha_y beta_y psi_y D_x Dprime_x D_y Dprime_y"
    for lf in lfinfo:
        print >>lffo, "%16s %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g"%(lf['name'], lf['s'], lf['alpha_x'], lf['beta_x'], lf['psi_x'], lf['alpha_y'],
                                                                                    lf['beta_y'], lf['psi_y'], lf['D_x'], lf['Dprime_x'], lf['D_y'], lf['Dprime_y'])

    lffo.close()

def do_lfcompare(lfinfo, lattice):
    lfinfo_len = len(lfinfo)
    lattice_len = len(lattice.get_elements())
    # first check number entries
    if lfinfo_len != lattice_len:
        print "# entries does not match."
        print "twiss file: ", lfinfo_len, " entries != Synergia/CHEF: ", lattice_len
        return

    # make lattice functions.  copied out of lattice_fns.py
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)
    syn_lf = []
    for elem in lattice.get_elements():
        syn_lf.append(lattice_simulator.get_lattice_functions(elem))

    i = 0
    for testfns in zip(lfinfo, syn_lf):
        madlf = testfns[0]
        syncheflf = testfns[1]
        if not test_close(madlf['s'], syncheflf.arc_length):
            print "Element ", i, " arclength MAD8: ", madlf['s'], " Synergia/CHEF: ", syncheflf.arc_length
        if not test_close(madlf['alpha_x'], syncheflf.alpha_x):
            print "Element ", i, " alpha_x MAD8: ", madlf['alpha_x'], " Synergia/CHEF: ", syncheflf.alpha_x
        if not test_close(madlf['alpha_x'], syncheflf.alpha_x):
            print "Element ", i, " alpha_x MAD8: ", madlf['alpha_x'], " Synergia/CHEF: ", syncheflf.alpha_x

        if not test_close(madlf['beta_x'], syncheflf.beta_x):
            print "Element ", i, " beta_x MAD8: ", madlf['beta_x'], " Synergia/CHEF: ", syncheflf.beta_x

        if not test_close(madlf['psi_x'], syncheflf.psi_x):
            print "Element ", i, " psi_x MAD8: ", madlf['psi_x'], " Synergia/CHEF: ", syncheflf.psi_x

        if not test_close(madlf['alpha_y'], syncheflf.alpha_y):
            print "Element ", i, " alpha_y MAD8: ", madlf['alpha_y'], " Synergia/CHEF: ", syncheflf.alpha_y

        if not test_close(madlf['beta_y'], syncheflf.beta_y):
            print "Element ", i, " beta_y MAD8: ", madlf['beta_y'], " Synergia/CHEF: ", syncheflf.beta_y

        if not test_close(madlf['psi_y'], syncheflf.psi_y):
            print "Element ", i, " psi_y MAD8: ", madlf['psi_y'], " Synergia/CHEF: ", syncheflf.psi_y

        if not test_close(madlf['D_x'], syncheflf.D_x):
            print "Element ", i, " D_x MAD8: ", madlf['D_x'], " Synergia/CHEF: ", syncheflf.D_x

        if not test_close(madlf['Dprime_x'], syncheflf.Dprime_x):
            print "Element ", i, " Dprime_x MAD8: ", madlf['Dprime_x'], " Synergia/CHEF: ", syncheflf.Dprime_x

        if not test_close(madlf['D_y'], syncheflf.D_y):
            print "Element ", i, " D_y MAD8: ", madlf['D_y'], " Synergia/CHEF: ", syncheflf.D_y

        if not test_close(madlf['Dprime_y'], syncheflf.Dprime_y):
            print "Element ", i, " Dprime_y MAD8: ", madlf['Dprime_y'], " Synergia/CHEF: ", syncheflf.Dprime_y


if __name__ == "__main__":
    try:
        option_values, leftover = getopt.getopt(sys.argv[1:], '', [
            'xmlfile=', 'lfcsvfile=', 'lfnpfile=', 'help'])
    except getopt.GetoptError as err:
        print err
        usage()
        sys.exit(10)

    if "--help" in [ov[0] for ov in option_values]:
        usage()
        sys.exit(0)

    xmlfile = None
    lfcsvfile = None
    lfnpfile = None
    lfcompare = False
    for ov in option_values:
        if ov[0] == "--xmlfile":
            xmlfile = ov[1]
        elif ov[0] == "--lfcsvfile":
            lfcsvfile = ov[1]
        elif ov[0] == "--lfnpfile":
            lfnpfile = ov[1]
        elif ov[0] == "--lfcompare":
            lfcompare = True

    if not leftover or len(leftover)!=1:
        print "error, need TWISS file"
        usage()
        sys.exit(10)

    twissfo = open(leftover[0])
    lfinfo, lattice = read_twiss_file(twissfo)
    twissfo.close()

    print "read lattice, number of elements: ", len(lattice.get_elements())
    print "lattice length: ", lattice.get_length()

    if xmlfile:
        synergia.lattice.xml_save_lattice(lattice, xmlfile)

    if lfnpfile:
        print "saving lattice functions to ", lfnpfile
        np.save(lfnpfile, lfinfo)

    if lfcsvfile:
        print "saving csv lattice functions to ", lfcsvfile
        write_lfcsvfile(lfcsvfile, lfinfo)

    if lfcompare:
        do_lfcompare(lfinfo, lattice)

    #print "Parsed lattice"
    #lattice.print_()

    sys.exit(0)
