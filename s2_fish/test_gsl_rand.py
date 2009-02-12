#/usr/bin/env python

import numpy
import gsl_rand
import time
import numpy

gr = gsl_rand.GSL_rand()

times = []
for i in range(0,10):
    t0 = time.time()
    a = gr.get_array((6,100000))
    t = time.time()
    times.append(t-t0)
print "time =",numpy.mean(times),"+/-",numpy.std(times)
