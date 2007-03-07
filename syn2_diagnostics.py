#!/usr/bin/env python
import Numeric
import math
import sys

x = 0
xprime = 1
y = 2
yprime = 3
z = 4
zprime = 5

x_xprime = 0
x_y = 1
xprime_y = 2
x_yprime = 3
xprime_yprime = 4
y_yprime = 5
x_z = 6
xprime_z = 7
y_z = 8
yprime_z = 9
x_zprime = 10
xprime_zprime = 11
y_zprime = 12
yprime_zprime = 13
z_zprime = 14

xprime_x = 0
y_x = 1
y_xprime = 2
yprime_x = 3
yprime_xprime = 4
yprime_y = 5
z_x = 6
z_xprime = 7
z_y = 8
z_yprime = 9
zprime_x = 10
zprime_xprime = 11
zprime_y = 12
zprime_yprime = 13
zprime_z = 14

def get_diagnostics(bunch):
    means = Numeric.zeros([6],Numeric.Float)
    stds = Numeric.zeros([6],Numeric.Float)
    for i in range(0,6):
        means[i] = Numeric.average(bunch[i,:])
    mom2s = Numeric.zeros([6,6],Numeric.Float)
    corrs = Numeric.zeros([6,6],Numeric.Float)
    for i in range(0,6):
        for j in range(i,6):
            tmp = bunch[i,:]*bunch[j,:]
            mom2s[i,j] = Numeric.average(tmp) - means[i]*means[j]
            mom2s[j,i] = mom2s[i,j]
    for i in range(0,6):
        stds[i] = math.sqrt(mom2s[i,i])
    for i in range(0,6):
        for j in range(i,6):
            corrs[i,j] = mom2s[i,j]/(stds[i]*stds[j])
            corrs[j,i] = corrs[i,j]
    return means,stds,mom2s,corrs

class Diagnostics:
    def __init__(self,units):
        self.s = []
        self.mean = [[],[],[],[],[],[]]
        self.std = [[],[],[],[],[],[]]
        self.u = units
        # n.b. we are ignoring mom2, corrs!

    def add(self,s,bunch):
        means,stds,mom2s,corrs = get_diagnostics(bunch.get_store().get_local_particles())
        self.s.append(s)
        for i in range(0,6):
            self.mean[i].append(means[i]/self.u[i])
            self.std[i].append(stds[i]/self.u[i])

    def get_coord_stds(self,bunch):
        means,stds,mom2s,corrs = get_diagnostics(bunch.get_store().get_local_particles())
        return (stds[x],stds[y],stds[z])
    
    def write(self,filename_prefix):
        # same format as fort.24
        f = open(filename_prefix+"_x.dat","w")
        for i in range(0,len(self.s)):
            f.write("%g %g %g %g %g\n" % \
                    (self.s[i],
                     self.mean[xprime][i], self.std[xprime][i],
                     self.mean[x][i], self.std[x][i]))
        f.close()
