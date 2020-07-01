#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

import tables

# get list of files
file_list = glob("turn_particles_0_*.h5")
file_list.sort()

# for f in file_list:
#     print "reading file:", f,
#     h5 = tables.open_file(f)
#     xstd = h5.root.particles[:, 0].std()
#     ystd = h5.root.particles[:, 2].std()
#     print " xstd: ", xstd, ", ystd: ", ystd
#     h5.close()

h5 = tables.open_file(file_list[0])
xstd = h5.root.particles[:, 0].std()
ystd = h5.root.particles[:, 2].std()
print "initial RMS xstd: ", xstd, ", ystd: ", ystd
h5.close()

h5diag = tables.open_file("full2_0.h5")
emitx = h5diag.root.emitx.read()
emity = h5diag.root.emity.read()
h5diag.close()

print "read diagnostics"

init_radius = 0.00415
cut4 = 4.0*init_radius
cut6 = 6.0*init_radius
print "cut4: ", cut4
print "cut4**2: ", cut4**2
print "cut6: ", cut6
print "cut6**2: ", cut6**2

cut4_particles = set()
cut6_particles = set()

ncut4 = [] # integrated # particles cut at 4 sigma
ncut6 = [] # integrated # particles cut at 6 sigma

xstdlist = []
xemitlist = []
xemitlist2 = []
turnlist = []

turn = 0
for f in file_list:
    print f, 
    h5 = tables.open_file(f)
    cut4list = (h5.root.particles[:,0]**2 + h5.root.particles[:,2]**2) >= cut4**2
    cut6list = (h5.root.particles[:,0]**2 + h5.root.particles[:,2]**2) >= cut6**2

    xcovar = np.cov(h5.root.particles[:,0], h5.root.particles[:,1])
    xemit = np.sqrt(xcovar[0,0]*xcovar[1,1] - xcovar[0,1]**2)
    xstdlist.append(h5.root.particles[:,0].std())
    xemitlist.append(xemit)
    turnlist.append(turn)

    add_to_list = [h5.root.particles[k, 6] for k in range(len(cut4list)) if cut4list[k]]
    cut4_particles.update(add_to_list)

    add_to_list = [h5.root.particles[k, 6] for k in range(len(cut6list)) if cut6list[k]]
    cut6_particles.update(add_to_list)

    ncut4.append(len(cut4_particles))
    ncut6.append(len(cut6_particles))
    print "turn: ", turn, "ncut4: ", ncut4[-1], " ncut6: ", ncut6[-1]

    turn += 25
    h5.close()

turns = np.array(turnlist)
print "end file loop"

np.savetxt("cut_particles.txt",
           np.vstack((turns, np.array(ncut4), np.array(ncut6))).transpose(),
                      header="# turn ncut4sigma ncut6sigma")

plt.figure()
plt.title("integrated particles cut at 4 sigma")
plt.plot(turns, np.array(ncut4))
plt.xlabel("turns")
plt.ylabel("particles cut")
plt.savefig("particles_cut_4sig.png")

plt.figure()
plt.title("integrated particles cut at 6 sigma")
plt.plot(turns, np.array(ncut6))
plt.xlabel("turns")
plt.ylabel("particles cut")
plt.savefig("particles_cut_6sig.png")

plt.figure()
plt.title("x std")
plt.plot(turns,np.array(xstdlist))
plt.xlabel("turns")
plt.ylabel("x std")
plt.savefig("x_std.png")

plt.figure()
plt.title("x RMS emittance")
plt.plot(turns, np.array(xemitlist))
plt.xlabel("turns")
plt.ylabel("x emit")
plt.savefig("x_emit.png")

h5_200 = tables.open_file("turn_particles_0_0008.h5")
r2_200 = h5_200.root.particles[:,0]**2 + h5_200.root.particles[:,2]**2
h5_200.close()
h5_800 = tables.open_file("turn_particles_0_0032.h5")
r2_800 = h5_800.root.particles[:,0]**2 + h5_800.root.particles[:,2]**2
h5_800.close()

plt.figure()
plt.title("r**2 for turn 200 and turn 800")
plt.hist(r2_800, bins=50, log=True, label='turn 800')
plt.hist(r2_200, bins=50, log=True, label='turn 200')
plt.xlabel("r**2")
plt.legend(loc='best')

plt.show()

    
