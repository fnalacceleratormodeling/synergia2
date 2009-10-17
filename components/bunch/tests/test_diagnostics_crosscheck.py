#!/usr/bin/env python

import numpy

total_num = 9;
particles = numpy.zeros([total_num,6],'d')
for part in range(0,total_num):
    for i in range(0,6):
        particles[part,i] = 10.0 * part + 1.1 * i
        
print "// python output for get_mean"
for i in range(0,6):
    print "BOOST_CHECK_CLOSE(diagnostics.get_mean()[%d], %0.15f, tolerance);" % \
        (i, numpy.mean(particles[:,i]))
print "// end python output for get_mean"
print

print "// python output for get_std"
for i in range(0,6):
    print "BOOST_CHECK_CLOSE(diagnostics.get_std()[%d], %0.15f, tolerance);" % \
        (i, numpy.std(particles[:,i]))
print "// end python output for get_std"
print



