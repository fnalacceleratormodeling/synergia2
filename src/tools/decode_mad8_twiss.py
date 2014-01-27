#!/usr/bin/env python

import sys
import os
from fortranformat import FortranRecordReader as frr
import getopt
import synergia

DEBUG = False

def usage():
    print "usage: decode_madx_twiss [options] TWISS-file"
    print "    options: --xmlfile=file write xml description of lattice"
    print "             --lfcsvfile=file write lattice functions in csv format to file"
    print "             --lfnpfile=file write numpy save file of lattice functions"
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
    for ov in option_values:
        if ov[0] == "--xmlfile":
            xmlfile = ov[1]
        elif ov[0] == "--lfcsvfile":
            lfcsvfile = ov[1]
        elif ov[0] == "lfnpfile":
            lfnpfile = ov[1]

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

    #print "Parsed lattice"
    #lattice.print_()

    sys.exit(0)
