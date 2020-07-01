#!/usr/bin/env python

import sys, os
import numpy as np
import scipy
import scipy.stats
from glob import glob
import tables
import matplotlib.pyplot as plt

file_list = glob("turn_particles_0_*.h5")
file_list.sort()

turns = []
kurt = []
minxs = []
maxxs = []

for f in file_list:
    h5 = tables.open_file(f)
    turn = h5.root.rep[()]
    if not f.endswith("_0_0000.h5"):
        turn += 25
    turns.append(turn)
    minxs.append(h5.root.particles[:,0].min())
    maxxs.append(h5.root.particles[:,0].max())
    kurt.append(scipy.stats.kurtosis(h5.root.particles[:,0]))
    h5.close()

plt.plot(turns, minxs, label="minimum x")
plt.plot(turns, maxxs, label="maximum x")
plt.xlabel("turn")
plt.ylabel("max or min")
plt.legend(loc='best')

plt.figure()
plt.plot(turns, kurt, label="x kurtosis")
plt.xlabel("turn")
plt.ylabel("kurtosis")
plt.legend(loc='best')

h5part0 = tables.open_file("turn_particles_0_0000.h5")
h5part100 = tables.open_file("turn_particles_0_0004.h5")
h5part1000 = tables.open_file("turn_particles_0_0040.h5")

kurt0 = scipy.stats.kurtosis(h5part0.root.particles[:,0])
kurt100 = scipy.stats.kurtosis(h5part100.root.particles[:,0])
kurt1000 = scipy.stats.kurtosis(h5part1000.root.particles[:,0])

plt.figure()
plt.hist(h5part100.root.particles[:,0], 100, [-0.022, 0.022], alpha=0.5, label="turn 100, kurtosis = %6.4f"%kurt100)
plt.hist(h5part0.root.particles[:,0], 100, [-0.022, 0.022], alpha=0.5, label="turn 0, kurtosis = %6.4f"%kurt0)
plt.xlabel("x")
plt.ylabel("density")
plt.legend(loc='best')

plt.figure()
plt.hist(h5part1000.root.particles[:,0], 100, [-0.022, 0.022], alpha=0.5, label="turn 1000, kurtosis = %6.4f"%kurt1000)
plt.hist(h5part0.root.particles[:,0], 100, [-0.022, 0.022], alpha=0.5, label="turn 0, kurtosis = %6.4f"%kurt0)
plt.xlabel("x")
plt.ylabel("density")
plt.legend(loc='best')


plt.show()

                
    
