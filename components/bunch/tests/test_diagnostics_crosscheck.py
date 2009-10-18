#!/usr/bin/env python

import numpy

total_num = 9;
particles = numpy.zeros([total_num, 6], 'd')
for part in range(0, total_num):
    for i in range(0, 6):
        particles[part, i] = 10.0 * part + (1.0 + part * part / 1000.0) * i
        
f = open("test_diagnostics_get_mean.icc", "w")
for i in range(0, 6):
    f.write("BOOST_CHECK_CLOSE(diagnostics.get_mean()[%d], %0.15f, tolerance);\n" % \
        (i, numpy.mean(particles[:, i])))
f.close()

f = open("test_diagnostics_get_std.icc", "w")
for i in range(0, 6):
    f.write("BOOST_CHECK_CLOSE(diagnostics.get_std()[%d], %0.15f, tolerance);\n" % \
        (i, numpy.std(particles[:, i])))
f.close()

mean = numpy.zeros([6], 'd')
sum2 = numpy.zeros([6, 6], 'd')
for i in range(0, 6):
    mean[i] = numpy.mean(particles[:, i])
for part in range(0, total_num):
    for i in range(0, 6):
        for j in range(0, 6):
            sum2[i, j] += (particles[part, i] - mean[i]) * \
                (particles[part, j] - mean[j])
mom2 = sum2 / total_num

f = open("test_diagnostics_get_mom2.icc", "w")
for i in range(0, 6):
    for j in range(0, 6):
        f.write("BOOST_CHECK_CLOSE(diagnostics.get_mom2()[%d][%d], %0.15f, tolerance);\n" % \
            (i, j, mom2[i, j]))
f.close()
