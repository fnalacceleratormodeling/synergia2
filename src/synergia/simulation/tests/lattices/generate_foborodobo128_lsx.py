#!/usr/bin/env python
import synergia
import numpy as np

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

lattice = synergia.lattice.Mad8_reader().get_lattice("model", "foborodobo128.lat")

# not setting the RF frequency, because that's what I'm testing
synergia.utils.write_lsexpr_file(lattice.as_lsexpr(), "foborodobo128_lattice.lsx")
