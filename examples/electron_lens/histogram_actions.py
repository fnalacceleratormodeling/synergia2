#!/usr/bin/env python

import sys, os
import numpy as np
import tables
import matplotlib.pyplot as plt

bpm_beta_x = 4.70792560035
bpm_alpha_x = -0.449832076233
bpm_beta_y = 46.1594967296
bpm_alpha_y = 3.37300523129

beta_x = 17.3259852015
alpha_x = 1.85063532729
beta_y = 17.2246414528
alpha_y = -1.90226103118

if len(sys.argv) < 2:
    print "usage: ", sys.argv[0], " file [file ...]"
    sys.exit(10)

os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

files = sys.argv[1:]
f1 = plt.figure()
plt.title("X actions")
f2 = plt.figure()
plt.title("Y actions")
xactions = []
yactions = []

for f in files:
    h5 = tables.open_file(f)
    x = h5.root.particles[:,0]
    xp = h5.root.particles[:,0]*alpha_x + \
               h5.root.particles[:,1]*beta_x
    xactions.append((x**2 + xp**2)/(2*beta_x))
    y = h5.root.particles[:, 2]
    yp = h5.root.particles[:,2]*alpha_y + \
               h5.root.particles[:,3]*beta_y
    yactions.append((y**2 + yp**2)/(2*beta_y))

plt.figure(1)
#plt.hist(xactions, bins=30, alpha=0.6, label=files)
for (xa, lbl) in zip(xactions, files):
    plt.hist(xa, alpha=0.5, bins=30,range=[0.0, 0.000025], label=lbl,log=True)
    #plt.hist(xa, alpha=0.5, bins=30, label=lbl,log=True)
plt.legend(loc='upper right')

plt.savefig("histogram_xactions.png")
plt.savefig("histogram_xactions.svg")

plt.figure(2)
for (ya, lbl) in zip(yactions, files):
    plt.hist(ya, alpha=0.5, bins=30, range=[0.0, 0.000025], label=lbl, log=True)
    #plt.hist(ya, alpha=0.5, bins=30, label=lbl, log=True)
#plt.hist(yactions, bins=30, label=files, alpha=0.6, range=[0.0, 0.20])
plt.legend(loc='upper right')

plt.savefig("histogram_yactions.png")
plt.savefig("histogram_yactions.svg")
plt.show()


