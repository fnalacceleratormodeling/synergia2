#!/usr/bin/env python3

import sys, os
import re
import numpy as np
import synergia

def remove_comments_and_continuation(lines):
    i = 0
    while i < len(lines):

        if re.search('!', lines[i]):
            # has a comment character, remove it and everything afterwards
            re.sub('!.*$', '\n', lines[i])

            # if it is totally whitespace now, delete the line
            if re.match('\s*', lines[i]):
                del lines[i]
            
            continue # keep going, I've dealt with the comment

        # check for a continuation, which could span multiple lines
        while len(lines[i]) > 2 and lines[i][-2] == '&':
            # add the next line on to the end of this one removing the comment character
            lines[i] = lines[i][:-2]+lines[i+1]
            del lines[i+1]
            
        # we can go on now
        i = i+1
        continue


# fixup problematic definitions in the MI lattice
# mi_lattice is a list of strings, each is one line in the lattice file
def mi_fixup(s):

    for i in range(len(s)):

        l = s[i]
        # LAM306 is a zero angle rbend.  replace with drift
        if l[0:7] == "LAM306:":
            l = 'LAM306: DRIFT, L=LILA_MAG\n'

        # K308[ABC] or 0 length rbends. replace with drift, with correct length
        if l[0:6] == "K308A:":
            l = "K308A: DRIFT, L=LRKB25_MAG\n"
        if l[0:6] == "K308B:":
            l = "K308B: DRIFT, L=LRKB25_MAG\n"
        if l[0:6] == "K308C:":
            l = "K308C: DRIFT, L=LRKD25_MAG\n" # different that K308A,B

        # fix the IDA RBEND magnets
        if re.search('(IDA\d\d\d:\s+)(RBEND, TYPE=IDA, L=LIDA_MAG, ANGLE=CIRBEND)', l):
            l = re.sub('L=LIDA_MAG', 'L=LIDA_MAG*(CIRBEND/2)/SIN(CIRBEND/2),E1=CIRBEND/2,E2=CIRBEND/2', l)
            l = re.sub(': RBEND', ': SBEND', l)
        # fix the IDB RBEND magnets
        if re.search('(IDB\d\d\d:\s+)(RBEND, TYPE=IDB, L=LIDB_MAG, ANGLE=CIRBEND)', l):
            l = re.sub('L=LIDB_MAG', 'L=LIDB_MAG*(CIRBEND/2)/SIN(CIRBEND/2),E1=CIRBEND/2,E2=CIRBEND/2', l)
            l = re.sub(': RBEND', ': SBEND', l)

        # IDC and IDD magnets have bending angle that is (2/3)*CIRBEND.
        # Notice that (1/2)(2/3)CIRBEND=CIRBEND/3.
        # fix the IDC RBEND magnets
        if re.search('(IDC\d\d\d:\s+)(RBEND, TYPE=IDC, L=LIDC_MAG)', l):
            l = re.sub('L=LIDC_MAG', 'L=LIDC_MAG*(CIRBEND/3)/SIN(CIRBEND/3),E1=CIRBEND/3,E2=CIRBEND/3', l)
            l = re.sub(': RBEND', ': SBEND', l)
            l = re.sub('ANGLE=0.666666666667', 'ANGLE=(2.0/3.0)', l)
        # fix the IDD RBEND magnets
        if re.search('(IDD\d\d\d:\s+)(RBEND, TYPE=IDD, L=LIDD_MAG)', l):
            l = re.sub('L=LIDD_MAG', 'L=LIDD_MAG*(CIRBEND/3)/SIN(CIRBEND/3),E1=CIRBEND/3,E2=CIRBEND/3', l)
            l = re.sub(': RBEND', ': SBEND', l)
            l = re.sub('ANGLE=0.666666666667', 'ANGLE=(2.0/3.0)', l)

        s[i] = l

if __name__ == "__main__":
    f = open('mi20_ra_08182020.lat')
    #f = open('testlattice.lat')
    lattice_lines = f.readlines()
    print("read ", len(lattice_lines), ' lines')
    remove_comments_and_continuation(lattice_lines)
    g = open('foo.lat', 'w')
    g.writelines(lattice_lines)
    g.close()
    print(len(lattice_lines), ' lines after removing comments and continuation')
    #m8_reader = synergia.lattice.Mad8_reader()
    #m8_reader.parse_string('\n'.join(lattice_lines))

    mi_fixup(lattice_lines)
    
    g = open('bar.lat', 'w')
    g.writelines(lattice_lines)
    g.close()

    m8_reader = synergia.lattice.Mad8_reader()
    m8_reader.parse_string('\n'.join(lattice_lines))

    lattice = m8_reader.get_lattice('ring_p_q100')

    print("lattice # elements: ", len(lattice.get_elements()))
    print("lattice length: ", lattice.get_length())

    lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)

    mi_beamline = lattice_simulator.get_chef_lattice().get_beamline()
    f = open('mi_beamline.txt', 'w')
    print(synergia.lattice.chef_beamline_as_string(mi_beamline), file=f)
    f.close()
    
    (nux, nuy) = lattice_simulator.get_both_tunes()
    chrx = lattice_simulator.get_horizontal_chromaticity()
    chry = lattice_simulator.get_vertical_chromaticity()

    print('Transverse tunes: qx: {}, qy: {}'.format(nux, nuy))
    print('Chromaticities: dqx: {}, dqy: {}'.format(chrx, chry))

    synergia.utils.write_lsexpr_file(lattice.as_lsexpr(), "mi20_ra_08182020.lsx")

    s = 0.0
    lf = []
    for elem in lattice.get_elements():
        elemlf = lattice_simulator.get_lattice_functions(elem)
        a = {}
        a['s'] = elemlf.arc_length
        a['beta_x'] = elemlf.beta_x
        a['beta_y'] = elemlf.beta_y
        a['alpha_x'] = elemlf.alpha_x
        a['alpha_y'] = elemlf.alpha_y
        a['name'] = elem.get_name()
        a['length'] = elemlf.arc_length
        a['D_x'] = elemlf.D_x
        a['D_y'] = elemlf.D_y
        a['Dprime_x'] = elemlf.Dprime_x
        a['Dprime_y'] = elemlf.Dprime_y
        a['psi_x'] = elemlf.psi_x
        a['psi_y'] = elemlf.psi_y
        lf.append(a)

    a0 = dict(lf[-1])
    a0['s'] = 0.0
    a0['psi_x'] = 0.0
    a0['psi_y'] = 0.0

    milf = []
    milf.append(a0)
    
    milf = milf +  lf

    np.save("mi_lf.npy", milf)
