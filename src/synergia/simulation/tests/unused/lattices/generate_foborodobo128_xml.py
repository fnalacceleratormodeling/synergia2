#!/usr/bin/env python
import synergia
import numpy as np

#######################################################

lattice = synergia.lattice.Mad8_reader().get_lattice("model", "foborodobo128.lat")

# not setting rf frequency because that's what I'm testing

synergia.lattice.xml_save_lattice(lattice, "foborodobo128_lattice.xml")

