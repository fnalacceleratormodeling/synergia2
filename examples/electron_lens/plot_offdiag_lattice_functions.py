#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt

lf = np.loadtxt("offdiag_lattice_functions.csv", usecols=(1,2,3,4,5,6,7))
print lf.shape
# arc, bx, by, ax, ay, psix, psiy
plt.title("Model beta functions")
plt.plot(lf[:,0], lf[:,1], label=r'$\beta_x$', lw=2)
plt.plot(lf[:,0], lf[:,2], label=r'$\beta_y$', lw=2)
plt.xlabel('s')
plt.ylabel(r'$\beta$ [m]')
plt.legend(loc='best')
plt.savefig("offdiag_beta_functions.png")
plt.savefig("offdiag_beta_functions.svg")
plt.show()
