#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.pyplot as plt

#    element      arc[m]     beta_x[m]      beta_y[m]     alpha_x     alpha_y       psi_x      psi_y      k1

lfsyn = np.loadtxt("CS_lattice_functions.csv", usecols=(1,2,3,4,5,6,7,8))
# name s l betx alfx bety alfy

s = lfsyn[:,0]
beta_x = lfsyn[:,1]
beta_y = lfsyn[:,2]
alpha_x = lfsyn[:,3]
alpha_y = lfsyn[:,4]
psi_x = lfsyn[:,5]
psi_y = lfsyn[:,6]

plt.figure()
plt.title(r'$\beta_x$')
plt.plot(s, beta_x)
plt.xlabel('s')
plt.ylabel(r'$\beta$')

plt.figure()
plt.title(r'$\alpha_x$')
plt.plot(s, alpha_x)
plt.xlabel('s')
plt.ylabel(r'$\alpha$')

plt.figure()
plt.title(r'$\beta_y$')
plt.plot(s, beta_y)
plt.xlabel('s')
plt.ylabel(r'$\beta$')

plt.figure()
plt.title(r'$\alpha_y$')
plt.plot(s, alpha_y)
plt.xlabel('s')
plt.ylabel(r'$\alpha$')

# M_{22} matrix element is
#
#  \sqrt{\frac{{\beta_A}{\beta_B}}} ( -alpha_B \sin \psi + \cos \psi )
#
# lenses are at rows 3, 8, 15, ...
#

# start at the fourth fodo 

porig = 3 + 3*12
#porig = 8 + 3*12

psi_b = psi_x[porig]
beta_b = beta_x[porig]
alpha_b = alpha_x[porig]

idx = range(porig-12, porig+12)
n = len(idx)
fact1 = np.zeros(n, dtype='d')
fact2 = np.zeros(n, dtype='d')
psi = np.zeros(n, dtype='d')

for i,p in enumerate(idx):
    psi_a = psi_x[p]
    beta_a = beta_x[p]
    psidiff = psi_b - psi_a
    psi[i] = psidiff
    betafact = np.sqrt(beta_a/beta_b)
    fact1[i] = -betafact * alpha_b * np.sin(psidiff)
    fact2[i] = betafact * np.cos(psidiff)

plt.figure()
plt.title(r'$-\sqrt{\frac{\beta_B}{\beta_A}} \alpha_B \sin \psi $')
plt.plot(psi, fact1)
plt.xlabel(r'$\psi$')

plt.figure()
plt.title(r'$\sqrt{\frac{\beta_B}{\beta_A}} \cos \psi $')
plt.plot(psi, fact2)
plt.xlabel(r'$\psi$')

plt.figure()
#plt.title(r'$\sqrt{\frac{\beta_B}{\beta_A}}(-\alpha_b \sin \psi +  \cos \psi ) $')
plt.title(r'$M_{22}$ at lens location')
plt.plot(psi, fact1+fact2)
plt.xlabel(r'$\psi$')

plt.savefig("m22.png")

plt.show()

