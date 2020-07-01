#!/usr/bin/env python
import sys, os
import matplotlib.pyplot as plt

from emit_analyze import *

#print "looking at directory: ", sys.argv[1]

de = dir_emittances(sys.argv[1])
#print "dir(de): ", dir(de)
plt.plot(de.turns, de.xemitrms, label='xemit RMS')
plt.plot(de.turns, de.xemit999, label='xemit 99.9%')
plt.legend(loc='best')
plt.figure()

plt.figure()
plt.plot(de.turns, de.betas_x, label='beta x')
plt.legend(loc='best')

plt.show()

