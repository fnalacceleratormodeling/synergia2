#!/usr/bin/env python
import synergia
import numpy as np

covariance = np.zeros((6,6), dtype='d')

for i in range(6):
    covariance[i,i] = float(i)
    if i%2 == 0:
        covariance[i,i+1] = float(i+0.5)
    else:
        covariance[i, i-1] = float(i-0.5)

synergia.convertors.xml_save_array2d(covariance,"covariance.xml")

print(np.array2string(covariance))
