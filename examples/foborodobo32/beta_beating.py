#!/usr/bin/env python
import sys
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt

import synergia

from beta_beating_options import opts

#####################################

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        raise RuntimeError, "map is unstable"

    mu =np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

#######################################################
# get focussing and defocussing quadrupoles for adjust tune
def get_fd_quads(lattice):
    f_quads = []
    d_quads = []
    for elem in lattice.get_elements():
        if elem.get_type() == "quadrupole":
            if elem.get_double_attribute("k1") > 0.0:
                f_quads.append(elem)
            elif elem.get_double_attribute("k1") < 0.0:
                d_quads.append(elem)
    return (f_quads, d_quads)

################################################################################

#CS_lattice entries
#    element      arc[m]     beta_x[m]      beta_y[m]     alpha_x     alpha_y       psi_x      psi_y       D_x[m]      D_y[m]      Dprime_x     Dprime_y
#

class LF:
    def __init__(self, name, arc, beta_x, beta_y, alpha_x, alpha_y, psi_x, psi_y, D_x, D_y, Dp_x, Dp_y):
        self.name = name
        self.arc = arc
        self.beta_x = beta_x
        self.beta_y = beta_y
        self.alpha_x = alpha_x
        self.alpha_y = alpha_y
        self.psi_x = psi_x
        self.psi_y = psi_y
        self.D_x = D_x
        self.D_y = D_y
        self.Dp_x = Dp_x
        self.Dp_y = Dp_y

def read_lattice_functions(lf_file):
    lattice_functions = []
    lf = open(lf_file)
    # readlines, skip first 2
    lines = lf.readlines()[2:]
    lf.close()

    for l in lines:
        fields = l.split()
        name = fields[0]
        arc = float(fields[1])
        beta_x = float(fields[2])
        beta_y = float(fields[3])
        alpha_x = float(fields[4])
        alpha_y = float(fields[5])
        psi_x = float(fields[6])
        psi_y = float(fields[7])
        D_x = float(fields[8])
        D_y = float(fields[9])
        Dp_x = float(fields[10])
        Dp_y = float(fields[11])
        lattice_functions.append(LF(name, arc, beta_x, beta_y, alpha_x, alpha_y, psi_x, psi_y, D_x, D_y, Dp_x, Dp_y))
    return lattice_functions

################################################################################


################################################################################
logger = synergia.utils.Logger(0)

# read the lattice in from a MadX sequence file
lattice = synergia.lattice.MadX_reader().get_lattice("model", "foborodobo32.madx")

# set the momentum of the beam.  This could have been in a beam statement
# in the lattice file, but this gives the option of setting it on the
# command line by creating a reference particle with the desired momentum.

# create the Reference particle object with the correct momentum
momentum = opts.momentum
four_momentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.mp)
four_momentum.set_momentum(momentum)

refpart = synergia.foundation.Reference_particle(1, four_momentum)

# set it into the lattice object
lattice.set_reference_particle(refpart)

energy = refpart.get_total_energy()
momentum = refpart.get_momentum()
gamma = refpart.get_gamma()
beta = refpart.get_beta()

print >>logger, "energy: ", energy
print >>logger, "momentum: ", momentum
print >>logger, "gamma: ", gamma
print >>logger, "beta: ", beta

length = lattice.get_length()
print >>logger, "lattice length: ", length

# set RF parameters. The RF voltage and phase (lag) are set to give a
# synchrotron tune and a stable bucket.
harmon = 32.0
freq = harmon * beta * synergia.foundation.pconstants.c/length

rf_volt = opts.rf_volt
rf_lag = opts.lag

print >>logger, "RF frequency: ", freq, " Hz"
print >>logger, "RF voltage: ", rf_volt, " MV"

for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", rf_volt)
        elem.set_double_attribute("lag", rf_lag)
        elem.set_double_attribute("freq", freq*1.0e-6)
        print >>logger, "rfcavity: ", elem.print_()

f = synergia.utils.Logger(0, "foborodobo32_lattice.out")
print >>f, lattice.as_string()
del f

# the lattice_simulator object lets us do some computations for
# lattice functions and other parameters.
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)

print >>logger, "created lattice_simulator"

myrank = 0
map = lattice_simulator.get_linear_one_turn_map()
print >>logger, "one turn map from synergia2.5 infrastructure"
print >>logger, np.array2string(map, max_line_width=200)

[l, v] = np.linalg.eig(map)

#print "l: ", l
#print "v: ", v

print >>logger, "eigenvalues: "
for z in l:
    print >>logger, "|z|: ", abs(z), " z: ", z, " tune: ", np.log(z).imag/(2.0*np.pi)

[ax, bx, qx] = map2twiss(map[0:2,0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6,4:6])

print >>logger, "Lattice parameters (assuming uncoupled map)"
print >>logger, "alpha_x: ", ax, " alpha_y: ", ay
print >>logger, "beta_x: ", bx, " beta_y: ", by
print >>logger, "q_x: ", qx, " q_y: ", qy
print >>logger, "q_z: ", qz, " beta_z: ", bz

alpha_c = lattice_simulator.get_momentum_compaction()
slip_factor = alpha_c - 1/gamma**2
print >>logger, "alpha_c: ", alpha_c, ", slip_factor: ", slip_factor

f_quads, d_quads = get_fd_quads(lattice)
print "len(f_quads): ", len(f_quads), " len(d_quads): ", len(d_quads)

(orig_xtune, orig_ytune) = lattice_simulator.get_both_tunes()
print >>logger, "Original base tunes, x: ", orig_xtune, " y: ", orig_ytune

# adjust tunes if requested by the xtune and/or the ytune parameter using
# the list of focussing or defocussing quadruples as adjusters.

do_adjust_tunes = False
if opts.xtune or opts.ytune:
    do_adjust_tunes = True
    if opts.xtune:
        target_xtune = opts.xtune
    else:
        target_xtune = orig_xtune
    if opts.ytune:
        target_ytune = opts.ytune
    else:
        target_ytune = orig_ytune

if do_adjust_tunes:
    print >>logger, "adjusting tunes, x: ", opts.xtune," y: ", opts.ytune
    lattice_simulator.adjust_tunes(target_xtune, target_ytune, f_quads, d_quads, 1.0e-6, 1)

hchrom = lattice_simulator.get_horizontal_chromaticity()
vchrom = lattice_simulator.get_vertical_chromaticity()

print >>logger, "horizontal chromaticity: %.16g"%hchrom
print >>logger, "vertical chromaticity: %.16g"%vchrom

lattice_simulator.print_lattice_functions()

shutil.move("CS_lattice_functions.dat", "CS_lf_orig.dat")

# increase quad strength 5% for the quads in the first 2 FODO cells

fquad = 0
dquad = 0
for elem in lattice.get_elements():
    if elem.get_type() == "quadrupole":
        k1 = elem.get_double_attribute("k1")
        if k1 > 0.0:
            fquad = fquad + 1
            elem.set_double_attribute("k1", 1.05*k1)

        if k1 < 0.0:
            dquad = dquad + 1
            elem.set_double_attribute("k1", 1.05*k1)

        if (fquad == 2) and (dquad == 2):
            break

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)
lattice_simulator.print_lattice_functions()

shutil.move("CS_lattice_functions.dat", "CS_lf_pert.dat")

lf_orig = read_lattice_functions("CS_lf_orig.dat")
lf_pert = read_lattice_functions("CS_lf_pert.dat")

plt.figure()
plt.title("beta x")
plt.plot(np.array([lf.arc for lf in lf_orig]),
         np.array([lf.beta_x for lf in lf_orig]), label='orig')
plt.plot(np.array([lf.arc for lf in lf_pert]),
         np.array([lf.beta_x for lf in lf_pert]), label='pert')
plt.legend(loc='best')

plt.figure()
plt.title("beta y")
plt.plot(np.array([lf.arc for lf in lf_orig]),
         np.array([lf.beta_y for lf in lf_orig]), label='orig')
plt.plot(np.array([lf.arc for lf in lf_pert]),
         np.array([lf.beta_y for lf in lf_pert]), label='pert')
plt.legend(loc='best')

plt.figure()
plt.title("Dispersion x")
plt.plot(np.array([lf.arc for lf in lf_orig]),
         np.array([lf.D_x for lf in lf_orig]), label='orig')
plt.plot(np.array([lf.arc for lf in lf_pert]),
         np.array([lf.D_x for lf in lf_pert]), label='pert')
plt.legend(loc='best')

plt.figure()
plt.title("Dispersion y")
plt.plot(np.array([lf.arc for lf in lf_orig]),
         np.array([lf.D_y for lf in lf_orig]), label='orig')
plt.plot(np.array([lf.arc for lf in lf_pert]),
         np.array([lf.D_y for lf in lf_pert]), label='pert')
plt.legend(loc='best')

plt.show()


