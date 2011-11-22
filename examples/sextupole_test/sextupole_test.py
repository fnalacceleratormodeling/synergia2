#!/usr/bin/env python
import synergia
import numpy as np
import beamline
import sys
# read in lattice containing a single thin sextupole
lattice1 = synergia.lattice.Mad8_reader().get_lattice("mylat1","sextupole.lat")
refpart1 = lattice1.get_reference_particle()
momentum1 = refpart1.get_momentum()
brho1 = (1.0e9/synergia.foundation.pconstants.c) * momentum1
# locate the sextupole and get its strength
sstr1 = -999999.0
for elem in lattice1.get_elements():
    if elem.get_name() == "mysex":
        sstr1 = elem.get_double_attribute("k2")
        break
if sstr1 == -999999.0:
    raise RuntimeError, "error, couldn't find sextupole element in lattice"
print "thin sextupole strength: ", sstr1
print "Brho: ", brho1

# need at least order 2 for sextupoles
lattice_simulator1 = synergia.simulation.Lattice_simulator(lattice1, 3)

# Read through the lattice and extract the sextupole strength
chefstr1 = -999999.0
cbmln = lattice_simulator1.get_chef_lattice().get_beamline()
for elm in cbmln:
    if elm.Name() == "mysex":
        chefstr1 = elm.Strength()
        break

if chefstr1 == -999999.0:
    raise RuntimeError, "error, couldn't find sextupole element in lattice"
print "Chef sextupole element strength: ", chefstr1

# read in lattice containing a single multitupole with a sextupole component
lattice2 = synergia.lattice.Mad8_reader().get_lattice("mylat2","mpsextupole.lat")
# need at least order 2 for sextupoles
lattice_simulator2 = synergia.simulation.Lattice_simulator(lattice2, 3)

# coordinate values to
xcoords = np.array([-2.0e-4, -1.0e-4, 0.0, 1.0e-4, 2.0e-4],'d')
xpcoords = np.array([-2.0e-5, -1.0e-5, 0.0, 1.0e-5, 2.0e-5],'d')
ycoords = np.array([-2.0e-4, -1.0e-4, 0.0, 1.0e-4, 2.0e-4],'d')
ypcoords = np.array([-2.0e-5, -1.0e-5, 0.0, 1.0e-5, 2.0e-5],'d')
zcoords = np.array([0.0],'d')
zpcoords = np.array([0.0],'d')

npart_test_bunch = len(xcoords) * len(xpcoords) * len(ycoords) * len(ypcoords) * len(zcoords) * len(zpcoords)

# construct a bunch
bunch1 = synergia.bunch.Bunch(lattice1.get_reference_particle(), npart_test_bunch, 1.0e11, synergia.MPI.COMM_WORLD)
                    
particles1 = bunch1.get_local_particles()

bunch2 = synergia.bunch.Bunch(lattice2.get_reference_particle(), npart_test_bunch, 1.0e11, synergia.MPI.COMM_WORLD)
                    
particles2 = bunch2.get_local_particles()


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

                        particles2[pindx,0] = x
                        particles2[pindx,1] = xp
                        particles2[pindx,2] = y
                        particles2[pindx,3] = yp
                        particles2[pindx,4] = z
                        particles2[pindx,5] = zp

                        pindx = pindx + 1

particles_orig = particles1.copy()

stepper1 = synergia.simulation.Independent_stepper_elements(
    lattice_simulator1, 1)
stepper2 = synergia.simulation.Independent_stepper_elements(
    lattice_simulator2, 1)

propagator1 = synergia.simulation.Propagator(stepper1)
propagator2 = synergia.simulation.Propagator(stepper2)

b1stepmd = synergia.bunch.Multi_diagnostics()
b1turnmd = synergia.bunch.Multi_diagnostics()
b1turnmd.append(synergia.bunch.Diagnostics_full2(bunch1, "bunch1.h5"))
propagator1.propagate(bunch1, 1, b1stepmd, b1turnmd, True)

b2stepmd = synergia.bunch.Multi_diagnostics()
b2turnmd = synergia.bunch.Multi_diagnostics()
b2turnmd.append(synergia.bunch.Diagnostics_full2(bunch2, "bunch2.h5"))
propagator2.propagate(bunch2, 1, b2stepmd, b2turnmd, True)

# is the old particles1 still a pointer to the bunch1 particles?
for i in range(npart_test_bunch):
    if (particles1[i,0] != bunch1.get_local_particles()[i,0]) or \
       (particles1[i,1] != bunch1.get_local_particles()[i,1]) or \
       (particles1[i,2] != bunch1.get_local_particles()[i,2]) or \
       (particles1[i,3] != bunch1.get_local_particles()[i,3]) or \
       (particles1[i,4] != bunch1.get_local_particles()[i,4]) or \
       (particles1[i,5] != bunch1.get_local_particles()[i,5]) :
       print "Error, particles1 is no longer a pointer to bunch1.get_local_particles()!!!"

# Make sure sextupole is the same as multipole
for i in range(npart_test_bunch):
    if (bunch1.get_local_particles()[i,0] != bunch2.get_local_particles()[i,0]) or \
       (bunch1.get_local_particles()[i,1] != bunch2.get_local_particles()[i,1]) or \
       (bunch1.get_local_particles()[i,2] != bunch2.get_local_particles()[i,2]) or \
       (bunch1.get_local_particles()[i,3] != bunch2.get_local_particles()[i,3]) or \
       (bunch1.get_local_particles()[i,4] != bunch2.get_local_particles()[i,4]) or \
       (bunch1.get_local_particles()[i,5] != bunch2.get_local_particles()[i,5]) :
       print "Error, bunch1 particles and different than bunch2 particles!!!"

# check that the sextupole has the correct effect
kickerror = False
for i in range(npart_test_bunch):
    x = particles1[i,0]
    y = particles1[i,2]
    pxkick = -(chefstr1/brho1) * (x**2 - y**2)
    pykick = (chefstr1/brho1) * 2.0*x*y
    apparent_pxkick = particles1[i,1]-particles_orig[i,1]
    apparent_pykick = particles1[i,3]-particles_orig[i,3]
    if abs(apparent_pxkick - pxkick) > 1.0e-15:
        print "Error: pxkick incorrect for particle at (x,y) = (",x,y,")"
        print "    propagate: ", apparent_pxkick," k*(x**2-y**2): ", pxkick
        kickerror = True
    if abs(apparent_pykick - pykick) > 1.0e-15:
        print "Error: pykick incorrect for particle at (x,y) = (",x,y,")"
        print "    propagate: ", apparent_pxkick, " k*2*x*y: ", pykick
        kickerror = True

if kickerror:
    print "Error conditions found"
    sys.exit(30)
else:
    print "no errors"
    sys.exit(0)
