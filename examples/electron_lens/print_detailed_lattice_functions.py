#!/usr/bin/env python

import sys, os
import numpy as np
import synergia
import matplotlib.pyplot as plt

#lattice_file = "el_Xoforodo.lsx"
#lattice_file = "el_Xoffdiag_yuri.lat"
lattice_file = "el_idealcomp.lat"

#lattice = synergia.lattice.Lattice(synergia.utils.read_lsexpr_file(lattice_file))
lattice = synergia.lattice.Mad8_reader().get_lattice("model", lattice_file)

fquads = []
dquads = []
for elem in lattice.get_elements():
    if elem.get_type() == "quadrupole" and elem.get_double_attribute("k1") > 0.0:
        fquads.append(elem)
    if elem.get_type() == "quadrupole" and elem.get_double_attribute("k1") < 0.0:
        dquads.append(elem)

print "there are ", len(fquads), " focussing quads"
print "there are ", len(dquads), " defocussing quads"

stepper = synergia.simulation.Independent_stepper_elements(lattice, 1, 1)
lattice_simulator = stepper.get_lattice_simulator()

lattice_simulator.adjust_tunes(0.74, 0.74, fquads, dquads, 1.0e-6)
lattice_simulator.update()
for elem in lattice.get_elements():
    if elem.get_type() == "quadrupole":
         print "quadrupole, k1=", elem.get_double_attribute("k1")

print "lattice functions at elens elements"
for elem in lattice.get_elements():
    if elem.get_type() == "elens":
        lf = lattice_simulator.get_lattice_functions(elem)
        print "arc_length: ", lf.arc_length, ", beta_x: ", lf.beta_x, ", alpha_x: ", lf.alpha_x, ", beta_y: ", lf.beta_y, ", alpha_y: ", lf.alpha_y

print
print "lattice functions at beginning/end of line:"
lf = lattice_simulator.get_lattice_functions(lattice.get_elements()[-1])
print "arc_length: ", lf.arc_length, ", beta_x: ", lf.beta_x, ", alpha_x: ", lf.alpha_x, ", beta_y: ", lf.beta_y, ", alpha_y: ", lf.alpha_y


(nux, nuy) = lattice_simulator.get_both_tunes()
print "tunes: ", nux, nuy

lattice_simulator.print_lattice_functions()

stepper = synergia.simulation.Independent_stepper(lattice, 1, 288)
stepper.force_update_operations_no_collective()
#stepper.print_()

ss = []
beta_x = []
beta_y = []
alpha_x = []
alpha_y = []

for step in stepper.get_steps():
    for o in step.get_operators():
        #o.print_()
        #print dir(o)
        slices = o.get_slices()
        lf = stepper.get_lattice_simulator().get_lattice_functions(slices[-1])
        ss.append(lf.arc_length)
        beta_x.append(lf.beta_x)
        beta_y.append(lf.beta_y)
        alpha_x.append(lf.alpha_x)
        alpha_y.append(lf.alpha_y)

plt.figure()
plt.title("beta functions")
plt.plot(ss, beta_x, label='beta x')
plt.plot(ss, beta_y, label='beta y')
plt.ylabel("beta")
plt.xlabel("s")
plt.legend(loc='best')

plt.figure()
plt.title("alpha functions")
plt.plot(ss, alpha_x, label='alpha x')
plt.plot(ss, alpha_y, label='alpha y')
plt.xlabel('s')
plt.ylabel('alpha')
plt.legend(loc='best')
plt.show()

sys.exit(0)
for elem in lattice.get_elements():
    if elem.get_name() == "boc" or elem.get_name() == "mlens":
        lf = lattice_simulator.get_lattice_functions(elem)
        print elem.get_name(), lf.arc_length, lf.beta_x, lf.alpha_x, lf.beta_y, lf.alpha_y

sys.exit(0)
dummy = synergia.simulation.Dummy_collective_operator("foo")
#stepper = synergia.simulation.Independent_stepper_elements(lattice, 1, 1)
stepper = synergia.simulation.Split_operator_stepper_elements(lattice, 1, dummy, 1)
ls = stepper.get_lattice_simulator()

steps = stepper.get_steps()
for s in steps:
    ops = s.get_operators()
    for o in ops:
        #print o.get_name(), o.get_type()
        if o.get_type() == "independent":
            #print "operator: ", o.get_name()
            if o.get_name() != "first_half":
                continue
            slices = o.get_slices()
            for s in slices:
                if s.get_lattice_element().get_name() != "bpm":
                    continue
                print s.get_lattice_element().get_name(),
                lf = ls.get_lattice_functions(s)
                print "bpm lattice functions"
                print lf.arc_length, lf.beta_x, lf.alpha_x, lf.beta_y, lf.alpha_y

last_elem = lattice.get_elements()[-1]
lf = ls.get_lattice_functions(last_elem)
print "ring lattice functions"
print lf.arc_length, lf.beta_x, lf.alpha_x, lf.beta_y, lf.alpha_y
