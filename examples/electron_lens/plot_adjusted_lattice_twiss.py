#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pyplot as plt

lfmad = np.loadtxt("adjusted_lattice_twiss.txt", skiprows=47, usecols=(1,2,3,4,5,6))
lfsyn = np.loadtxt("adjusted_lattice_synergia_lf.txt", skiprows=47, usecols=(1,2,3,4,5,6))
# name s l betx alfx bety alfy

plt.figure()
plt.title("beta x comparison")
plt.plot(lfmad[:,0], lfmad[:,2], label='madx beta x')
plt.plot(lfsyn[:,0], lfsyn[:,2], label='synergia beta x')
plt.xlabel('s')
plt.ylabel('beta')
plt.legend(loc='best')

plt.figure()
plt.title("beta y comparison")
plt.plot(lfmad[:,0], lfmad[:,4], label='mad beta y')
plt.plot(lfsyn[:,0], lfsyn[:,4], label='synergia beta y')
plt.xlabel('s')
plt.ylabel('beta')
plt.legend(loc='best')

plt.figure()
plt.title("Model beta functions")
plt.plot(lfsyn[:,0], lfsyn[:,2], lw=2, label=r'$\beta_x$')
plt.plot(lfsyn[:,0], lfsyn[:,4], lw=2, label=r'$\beta_y$')
plt.xlabel('s [m]')
plt.ylabel(r'$\beta$ [m]')
plt.legend(loc='best', fontsize='x-large')
plt.savefig('adjusted_beta_functions.png')

plt.show()

