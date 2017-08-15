#!/usr/bin/env python
# -*- coding: utf-8 -*-
import synergia
from modes_options import opts

lattice = synergia.lattice.Mad8_reader().get_lattice("model", "modes_lattice.lat")
lattice_length=lattice.get_length() 
reference_particle = lattice.get_reference_particle()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
energy = reference_particle.get_total_energy()
harmon=opts.num_buckets
freq=harmon*beta*synergia.foundation.pconstants.c/lattice_length



#freq=37867099.7584
print "initial frequency=",freq, "Hz"    


for elem in lattice.get_elements():
    if opts.chef_propagate:
        elem.set_string_attribute("extractor_type", "chef_propagate")
    elif  opts.chef_map:   
        elem.set_string_attribute("extractor_type", "chef_map")


    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", opts.rf_voltage)
        elem.set_double_attribute("freq", freq)
        elem.set_double_attribute("lag", 0.)








synergia.lattice.xml_save_lattice(lattice, "modes_lattice.xml")

