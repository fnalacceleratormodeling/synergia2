#!/usr/bin/env python
import synergia
import numpy as np
import beamline
import sys

# read in lattice containing a single thin quadrupole
lattice1 = synergia.lattice.Mad8_reader().get_lattice("mylat1","quadrupole_with_sextupole.lat")
refpart1 = lattice1.get_reference_particle()
momentum1 = refpart1.get_momentum()
#print "dir(refpart1): ", dir(refpart1)
energy1 = refpart1.get_total_energy()

brho1 = (1.0e9/synergia.foundation.pconstants.c) * momentum1
# locate the quadrupole and get its strength

# propagage a particle through a drift using the full chef definition
# of coordinates.  particle is a list of [x, npx, y, npy, cdt, ndp]
def drift(particle, length):
    # get px, py, and pz
    ptot = momentum1 * (1.0 + particle[5])
    px = momentum1 * particle[1]
    py = momentum1 * particle[3]
    pz = np.sqrt(ptot**2 - px**2 - py**2)
    particle[0] += px/pz * length
    particle[2] += py/pz * length

lattice1.print_()

print "Brho: ", brho1
print "momentum: ", momentum1
print "energy: ", energy1

# need at least order 2 for sextupoles
lattice_simulator1 = synergia.simulation.Lattice_simulator(lattice1, 3)

#print "dir(lattice_simulator1): ", dir(lattice_simulator1)

#print "dir(lattice_simulator1.get_chef_lattice()): ", dir(lattice_simulator1.get_chef_lattice())

# need to create the stepper in order for the beam line to be sliced

stepper1 = synergia.simulation.Independent_stepper_elements(
    lattice_simulator1, 1)

sliced_beamline = lattice_simulator1.get_chef_lattice().get_sliced_beamline()

pr1 = beamline.Proton(energy1)

# find the quadrupole strength
quadstrengths = []
quadlengths = []
thinpole_strength = None

for elm in sliced_beamline:
    print elm.Name(), elm.Type(), elm.Strength(), elm.OrbitLength(pr1)

    if elm.Type() == "quadrupole":
        # I've found one of the quads
        quadstrengths.append(elm.Strength())
        quadlengths.append(elm.OrbitLength(pr1))
    elif elm.Type() == "ThinPole":
        thinpole_strength = elm.Strength()


if len(quadstrengths) == 0:
    raise RuntimeError, "error,didn't find exactly 2 quadrupole elements in lattice"
if not thinpole_strength:
    raise RuntimeError, "error, didn't find a ThinPole element in lattice"

print "CHEF quadrupole element strengths: ", quadstrengths
print "CHEF quadrupole lengths: ", quadlengths
print "CHEF ThinPole strength: ", thinpole_strength

# The thinpole strength should be equal the integrated B'*L of the quadrupoles

tp_shouldbe = np.vdot(quadlengths, quadstrengths)
if thinpole_strength != tp_shouldbe:
    raise RuntimeError, "error: thin_pole strength: %g incorrect, should be %g"%(thinpole_strength, tp_shouldbe)

# set the quad strengths to 0
for elm in sliced_beamline:
    if elm.Type() == "quadrupole":
        elm.setStrength(0.0)

# Now print out beamline again
print "after setting strengths to 0"
for elm in sliced_beamline:
    print elm.Name(), elm.Type(), elm.Strength(), elm.OrbitLength(pr1)

pr2 = beamline.Proton(energy1)
pr2.set_x(.00001)
pr2.set_npx(0.0)
pr2.set_y(-0.00002)
pr2.set_npy(0.0)

print "propagating particle: ", pr2.get_x(), pr2.get_npx(), pr2.get_y(), pr2.get_npy(), pr2.get_cdt(), pr2.get_ndp()
for elm in sliced_beamline:
    print elm.Name(), elm.Type(), elm.Strength(), elm.OrbitLength(pr1)
    elm.propagate(pr2)
    print pr2.get_x(), pr2.get_npx(), pr2.get_y(), pr2.get_npy(), pr2.get_cdt(), pr2.get_ndp()

# coordinate values to test
xcoords = np.array([-2.0e-4, -1.0e-4, 0.0, 1.0e-4, 2.0e-4],'d')
xpcoords = np.array([-2.0e-5, -1.0e-5, 0.0, 1.0e-5, 2.0e-5],'d')
ycoords = np.array([-2.0e-4, -1.0e-4, 0.0, 1.0e-4, 2.0e-4],'d')
ypcoords = np.array([-2.0e-5, -1.0e-5, 0.0, 1.0e-5, 2.0e-5],'d')
zcoords = np.array([0.0],'d')
zpcoords = np.array([0.0],'d')

npart_test_bunch = len(xcoords) * len(xpcoords) * len(ycoords) * len(ypcoords) * len(zcoords) * len(zpcoords)

# construct a bunch
bunch1 = synergia.bunch.Bunch(lattice1.get_reference_particle(), npart_test_bunch, 1.0e11, synergia.utils.Commxx())

particles1 = bunch1.get_local_particles()

pindx = 0
for x in xcoords:
    for xp in xpcoords:
        for y in ycoords:
            for yp in ypcoords:
                for z in zcoords:
                    for zp in zpcoords:
                        particles1[pindx,0] = x
                        particles1[pindx,1] = xp
                        particles1[pindx,2] = y
                        particles1[pindx,3] = yp
                        particles1[pindx,4] = z
                        particles1[pindx,5] = zp

                        pindx = pindx + 1

particles_orig = particles1.copy()

propagator1 = synergia.simulation.Propagator(stepper1)

b1stepmd = synergia.bunch.Multi_diagnostics()
b1turnmd = synergia.bunch.Multi_diagnostics()
b1turnmd.append(synergia.bunch.Diagnostics_full2(bunch1, "bunch1.h5"))
propagator1.propagate(bunch1, 1, b1stepmd, b1turnmd, True)

# is the old particles1 still a pointer to the bunch1 particles?
for i in range(npart_test_bunch):
    if (particles1[i,0] != bunch1.get_local_particles()[i,0]) or \
       (particles1[i,1] != bunch1.get_local_particles()[i,1]) or \
       (particles1[i,2] != bunch1.get_local_particles()[i,2]) or \
       (particles1[i,3] != bunch1.get_local_particles()[i,3]) or \
       (particles1[i,4] != bunch1.get_local_particles()[i,4]) or \
       (particles1[i,5] != bunch1.get_local_particles()[i,5]) :
       print "Error, particles1 is no longer a pointer to bunch1.get_local_particles()!!!"

# check that the propagation worked as expected
kickerror = False
for i in range(npart_test_bunch):
    # get copy original particle
    p_orig = particles_orig[i,:].copy()
    drift(p_orig, quadlengths[0])
    x = p_orig[0]
    y = p_orig[2]
    # this is a little off because we're not using the full definition
    # of px/p_pref
    # b2 value is 0.001
    xpkick = -0.001 * (thinpole_strength/brho1) * (x**2 - y**2)
    ypkick =  0.001 * (thinpole_strength/brho1) * 2.0 * x * y
    p_orig[1] += xpkick
    p_orig[3] += ypkick
    drift(p_orig, quadlengths[1])

    apparent_pxkick = particles1[i,1]-particles_orig[i,1]
    apparent_pykick = particles1[i,3]-particles_orig[i,3]
    if abs(apparent_pxkick - xpkick) > 1.0e-15:
        print "Error: xpkick incorrect for particle at (x,y) = (",x,y,")"
        print "    propagate: ", apparent_pxkick," k*(x**2-y**2): ", xpkick
        kickerror = True
    if abs(apparent_pykick - ypkick) > 1.0e-15:
        print "Error: pykick incorrect for particle at (x,y) = (",x,y,")"
        print "    propagate: ", apparent_pykick, " 2*k*x*y: ", ypkick
        kickerror = True

if kickerror:
    print "Errors detected"
    sys.exit(30)
else:
    print "no errors"
    sys.exit(0)
