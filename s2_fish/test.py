#!/usr/bin/env python

from s2_fish import *

f = Scalar_Field()
f.set_num_points(int3(16,16,16))
i3 = int3(1,2,3)
print "i3[2] =",i3.get(2)
