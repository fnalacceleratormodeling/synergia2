#!/usr/bin/env python3
from __future__ import print_function
import sys
import os
import numpy as np
import tables
import matplotlib.pyplot as plt


tracks_file = "tracks.h5"
nf_file = "nf.dat"

# read in track data
h5 = tables.open_file(tracks_file,'r')
track_coords = h5.root.track_coords.read()
h5.close()

# read in normal form data
nf = np.loadtxt(nf_file)

pnum = int(sys.argv[1])

# input is (# turns, number_particles*(x, xp, y, yp, z, zp))
nfdata = nf[:, 6*pnum:6*(pnum+1)]

xphase = np.arctan2(nfdata[:,1],nfdata[:,0])
yphase = np.arctan2(nfdata[:,3],nfdata[:,2])
xtune = (xphase[1:] - xphase[0:-1])/(2*np.pi)
negmask = xtune < 0.0
xtune = xtune + negmask * 1
ytune = (yphase[1:] - yphase[0:-1])/(2*np.pi)
negmask = ytune < 0.0
ytune = ytune + negmask * 1



plt.figure()
plt.title("x vs. xp particle %d"%pnum)
plt.plot(track_coords[:, pnum, 0], track_coords[:, pnum, 1], 'o')
plt.xlabel("x")
plt.ylabel("xp")

plt.figure()
plt.title("y vs. yp particle %d"%pnum)
plt.plot(track_coords[:, pnum, 2], track_coords[:, pnum, 3], 'o')
plt.xlabel("y")
plt.ylabel("yp")

plt.figure()
plt.title("Re vs. Im particle %d"%pnum)
plt.plot(nfdata[:,0], nfdata[:,1], '.', label="Re(a0) vs Im(a0)")
plt.plot(nfdata[:,2], nfdata[:,3], '.', label="Re(a1) vs Im(a1)")
plt.xlabel("Re(a)")
plt.ylabel("Im(a)")
plt.legend(loc='best')

plt.figure()
plt.title("sqrt(I) particle %d"%pnum)
plt.plot(np.sqrt(nfdata[:,0]**2 + nfdata[:,1]**2), label="sqrt(I0)")
plt.plot(np.sqrt(nfdata[:,2]**2 + nfdata[:,3]**2), label="sqrt(I1)")
plt.legend(loc='best')
plt.xlabel("turn")
plt.ylabel("sqrt(action)")

plt.figure()
plt.title("angle 1 vs. angle 0")
plt.plot(xphase, yphase, '.')
plt.xlabel("angle 0")
plt.ylabel("angle 1")

plt.subplot(211)
plt.plot(xtune)
plt.ylabel("x tune")
plt.subplot(212)
plt.plot(ytune)
plt.ylabel("y tune")

plt.show()
