#!/usr/bin/env python

import populate
import Numeric
import MLab
from ArrayPrinter import array2string

p = Numeric.zeros((1000,7),'d')
covs = Numeric.zeros((6,6),'d')
means = Numeric.zeros((6,),'d')

for i in range(0,6):
    means[i] = 0.01*i;
    for j in range(0,6):
        if i==j:
            covs[i,j] = i+1.0
        if (i+1 == j) and (i%2 == 0):
            covs[i,j] = 0.1*(j+1.0)
            covs[j,i] = 0.1*(j+1.0)

print "requested covs"
print covs
print "requested means"
print means

populate.populate_6d_gaussian(p,means,covs,0)

cp = MLab.cov(p)
print array2string(cp[0:6,0:6],suppress_small=1,precision=4)
mp = MLab.mean(p)
print mp

print "covariances after"
print covs
