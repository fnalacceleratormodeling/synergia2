#!/usr/bin/env python

import sys, os
import numpy as np
import synergia

lattice = synergia.lattice.Lattice(synergia.utils.read_lsexpr_file("adjusted_lattice.lsx"))

#print "lattice: ", lattice.as_string()

s = 0.0
quadcnt = 0
mrkcnt = 0
rfcnt = 0
mpcnt = 0
elenscnt = 0

for elem in lattice.get_elements():
    elem.set_double_attribute("at", s+elem.get_length()/2.0)
    if elem.get_type() == "quadrupole":
        elem.set_string_attribute("newname", "q%04d"%quadcnt)
        quadcnt += 1
    elif elem.get_type() == "drift":
        pass
    elif elem.get_type() == "marker":
        elem.set_string_attribute("newname", "mrk%04d"%mrkcnt)
        mrkcnt += 1
    elif elem.get_type() == "elens":
        elem.set_string_attribute("newname", "lens%04d"%elenscnt)
        elenscnt += 1
    elif elem.get_type() == "multipole":
        elem.set_string_attribute("newname", "mp%04d"%mpcnt)
        mpcnt += 1
    elif elem.get_type() == "rfcavity":
        elem.set_string_attribute("newname", "rfc%04d"%rfcnt)
        rfcnt += 1
    s += elem.get_length()

refpart = lattice.get_reference_particle()
print "beam, particle=proton, energy=%.17g;"%refpart.get_total_energy()

for elem in lattice.get_elements():
    if elem.get_type() == "quadrupole":
        print elem.get_string_attribute("newname"), ": quadrupole, ",
        print "k1=%.17g, "%elem.get_double_attribute("k1"),
        print "l=%.17g;"%elem.get_length()
    elif elem.get_type() == "marker" or elem.get_type() == "elens":
        print elem.get_string_attribute("newname"), ": marker;"
    elif elem.get_type() == "multipole":
        print elem.get_string_attribute("newname"), ": multipole, knl={0.0};"
    elif elem.get_type() == "rfcavity":
        print elem.get_string_attribute("newname"), ": rfcavity, ",
        print "l=%.17g, "%elem.get_double_attribute("l"),
        print "freq=%.17g, "% elem.get_double_attribute("freq"),
        print "lag=%.17g, "% elem.get_double_attribute("lag"),
        print "volt=%.17g, "% elem.get_double_attribute("volt"),
        print "harmon=%d;"% int(elem.get_double_attribute("harmon"))

print
print "machine: sequence, l=", lattice.get_length(), ";"
for elem in lattice.get_elements():
    if ((elem.get_type() == "quadrupole") or
        (elem.get_type() == "marker") or
        (elem.get_type() == "elens") or
        (elem.get_type() == "multipole") or
        (elem.get_type() == "rfcavity")):
        print elem.get_string_attribute("newname"), ", at=%.17g;"% elem.get_double_attribute("at")

print "endsequence;"

stepper = synergia.simulation.Independent_stepper_elements(lattice, 1, 1)
lattice_simulator = stepper.get_lattice_simulator()

f = open("adjusted_lattice_synergia_lf.txt", "w")
# write 46 lines of stuff
for l in range(46):
    print >>f, "#intro line %d"%l
print >>f, "#  name, s, l, betx alfx bety alfy"
s = 0.0
elf = lattice_simulator.get_lattice_functions(lattice.get_elements()[-1])
print >>f, "start ", 0, 0.0, elf.beta_x, elf.alpha_x, elf.beta_y, elf.alpha_y
for elem in lattice.get_elements():
    elf = lattice_simulator.get_lattice_functions(elem)
    print >>f, elem.get_name(), " ", elf.arc_length, elem.get_length(), elf.beta_x, elf.alpha_x, elf.beta_y, elf.alpha_y
    s += elem.get_length()
print >>f, "end ", elf.arc_length, 0.0, elf.beta_x, elf.alpha_x, elf.beta_y, elf.alpha_y
f.close()
