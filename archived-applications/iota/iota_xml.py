#!/usr/bin/env python
# -*- coding: utf-8 -*-
import synergia
#from modes_options import opts

lattices = {}
lattice_repo = './'
lattices['t1_1IO_66'] = lattice_repo + "lattice_1IO_center.madx" #centered t1 6.6 1IO lattice
lattices['t3_1IO_66'] = lattice_repo + "lattice_1IO_nll_center.madx" #centered t3 6.6 1IO lattice

USE_NLL = False #Flag for using the version of the lattice with the nonlinear element
if USE_NLL:
    name = 't3_1IO_66'
else:
    name = 't1_1IO_66'

lattice = synergia.lattice.MadX_reader().get_lattice("iota", lattices[name])  



#freq=37867099.7584
#print "initial frequency=",freq, "Hz"    







xmllattice_name="iota_"+name +".xml"

synergia.lattice.xml_save_lattice(lattice, xmllattice_name)
print "lattice xml file saved in ",xmllattice_name

