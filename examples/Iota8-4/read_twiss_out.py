#!/usr/bin/env python
import sys, os
import numpy as np

# first 52 rows are header
# row 53 has contents:
# * NAME                                  POS                      L                   BETX                   ALFX                    MUX                   BETY                   ALFY                    MUY                   RE11                   RE12                   RE13                   RE14                   RE15                   RE16                   RE21                   RE22                   RE23                   RE24                   RE25                   RE26                   RE31                   RE32                   RE33                   RE34                   RE35                   RE36                   RE41                   RE42                   RE43                   RE44                   RE45                   RE46                   RE51                   RE52                   RE53                   RE54                   RE55                   RE56                   RE61                   RE62                   RE63                   RE64                   RE65                   RE66                     RE 
re = np.loadtxt('twiss.out', skiprows=52, usecols=range(9, 45)).reshape((6,6))

print('re: ')
print(np.array2string(re, max_line_width=200))
