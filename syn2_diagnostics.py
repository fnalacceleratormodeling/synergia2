#!/usr/bin/env python
import Numeric
import LinearAlgebra
from math import sqrt
import sys

import s2_diagnostics

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

def get_spatial_means_stds(bunch):
    means = Numeric.zeros([3],Numeric.Float)
    stds = Numeric.zeros([3],Numeric.Float)
    s2_diagnostics.get_spatial_means_stds(bunch.get_store(),means,stds)
    return means,stds

def old_get_diagnostics(bunch):
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
        stds[i] = sqrt(mom2s[i,i])
    for i in range(0,6):
        for j in range(i,6):
            corrs[i,j] = mom2s[i,j]/(stds[i]*stds[j])
            corrs[j,i] = corrs[i,j]
    return means,stds,mom2s,corrs

def get_diagnostics(bunch,units):
    means = Numeric.zeros([6],Numeric.Float)
    mom2s = Numeric.zeros([6,6],Numeric.Float)
    corrs = Numeric.zeros([6,6],Numeric.Float)
    diagmom4s = Numeric.zeros([6],Numeric.Float)
    s2_diagnostics.get_moments_corrs(bunch.get_store(),units,means,mom2s,corrs,diagmom4s)
    return means,mom2s,corrs,diagmom4s

class Diagnostics:
    def __init__(self,units):
        self.s = []
        self.means = []
        self.mom2s = []
        self.corrs = []
        self.diagmom4s = []
        self.stds = []
        self.emitxs = []
        self.emitys = []
        self.emitzs = []
        self.emitxys = []
        #~ self.emitxzs = []
        self.emitxyzs = []
        self.u = units

    def add(self,s,bunch):
        means,mom2s,corrs,diagmom4s = get_diagnostics(bunch,self.u)
        self.s.append(s)
        self.means.append(means)
        self.mom2s.append(mom2s)
        self.corrs.append(corrs)
        self.diagmom4s.append(diagmom4s)
        #derived quantities
        self.stds.append(Numeric.sqrt(Numeric.diagonal(mom2s)))
        self.emitxs.append(sqrt(abs(LinearAlgebra.determinant(mom2s[0:2,0:2]))))
        self.emitys.append(sqrt(abs(LinearAlgebra.determinant(mom2s[2:4,2:4]))))
        self.emitzs.append(sqrt(abs(LinearAlgebra.determinant(mom2s[4:6,4:6]))))
        self.emitxys.append(sqrt(abs(LinearAlgebra.determinant(mom2s[0:4,0:4]))))
        #~ self.emitxzs.append(sqrt(LinearAlgebra.determinant(
        self.emitxyzs.append(sqrt(abs(LinearAlgebra.determinant(mom2s))))
    
    def get_s(self):
        return Numeric.array(self.s)
    def get_means(self):
        return Numeric.array(self.means)
    def get_stds(self):
        return Numeric.array(self.stds)
     
    def write(self,filename_prefix):
        return
        # same format as fort.24
        fx = open(filename_prefix+"_x.dat","w")
        fy = open(filename_prefix+"_y.dat","w")
        fz = open(filename_prefix+"_z.dat","w")
        for i in range(0,len(self.s)):
            fx.write("%g %g %g %g %g\n" % \
                    (self.s[i],
                     self.mean[xprime][i], self.std[xprime][i],
                     self.mean[x][i], self.std[x][i]))
            fy.write("%g %g %g %g %g\n" % \
                    (self.s[i],
                     self.mean[yprime][i], self.std[yprime][i],
                     self.mean[y][i], self.std[y][i]))
            fz.write("%g %g %g %g %g\n" % \
                    (self.s[i],
                     self.mean[zprime][i], self.std[zprime][i],
                     self.mean[z][i], self.std[z][i]))
        fx.close()
        fy.close()
        fz.close()
