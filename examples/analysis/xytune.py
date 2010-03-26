#!/usr/bin/env python

# plot the spectrum of the coherent dipole oscillation
# from a synergia run.

# set lettering to serif font
from matplotlib import rcParams
rcParams['font.family'] = 'serif'

# import stuff we need

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import tables

if len(sys.argv) == 1:
    print "give me a file name!"
    sys.exit(10)

if not os.access(sys.argv[1], os.R_OK):
    print "I can't read ", sys.argv[1]
    sys.exit(10)

# open the data file
h5file = tables.openFile(sys.argv[1])

# pick out the mean position of x and y
xmean = h5file.root.mean[:,0]
ymean = h5file.root.mean[:,2]

# get spectrum
xspect = np.abs(np.fft.fft(xmean))
yspect = np.abs(np.fft.fft(ymean))

# n is the number of turns in the file
n = xmean.shape[0]
df = 1.0/n
f = np.arange(n,dtype='d')*df

plt.plot(f, xspect, 'r', f, yspect, 'g')
plt.legend(('x spectrum', 'y spectrum'))
ax = plt.gca()
plt.setp(ax, xlabel='tune', ylabel='strength')
plt.show()
