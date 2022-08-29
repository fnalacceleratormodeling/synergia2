#!/usr/bin/env python
from math import sqrt
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

random_particles = numpy.random.lognormal(size=[total_num,6])
corr = numpy.corrcoef(random_particles.transpose())
f = open("test_diagnostics_get_random_particles.icc", "w")
for i in range(0, total_num):
    for j in range(0, 6):
        f.write("particles[%d][%d] = %0.15f;\n" % \
            (i, j, random_particles[i,j]))
f.close()
f = open("test_diagnostics_get_corr.icc", "w")
for i in range(0, 6):
    for j in range(0, 6):
        f.write("BOOST_CHECK_CLOSE(diagnostics.get_corr()[%d][%d], %0.15f, tolerance_corr);\n" % \
            (i, j, corr[i,j]))
f.close()

mean = numpy.zeros([6], 'd')
sum2 = numpy.zeros([6, 6], 'd')
for i in range(0, 6):
    mean[i] = numpy.mean(random_particles[:, i])
for part in range(0, total_num):
    for i in range(0, 6):
        for j in range(0, 6):
            sum2[i, j] += (random_particles[part, i] - mean[i]) * \
                (random_particles[part, j] - mean[j])
random_mom2 = sum2 / total_num
#f = open("test_diagnostics_get_random_mom2.icc", "w")
#for i in range(0, 6):
#    for j in range(0, 6):
#        f.write("BOOST_CHECK_CLOSE(diagnostics.get_mom2()[%d][%d], %0.15f, tolerance);\n" % \
#            (i, j, random_mom2[i, j]))
#f.close()
f = open("test_diagnostics_get_emitx.icc","w")
f.write("BOOST_CHECK_CLOSE(diagnostics.get_emitx(), %0.15f, tolerance_emit2d);\n" % \
       sqrt(numpy.linalg.det(random_mom2[0:2,0:2])))
f.close()
f = open("test_diagnostics_get_emity.icc","w")
f.write("BOOST_CHECK_CLOSE(diagnostics.get_emity(), %0.15f, tolerance_emit2d);\n" % \
       sqrt(numpy.linalg.det(random_mom2[2:4,2:4])))
f.close()
f = open("test_diagnostics_get_emitz.icc","w")
f.write("BOOST_CHECK_CLOSE(diagnostics.get_emitz(), %0.15f, tolerance_emit2d);\n" % \
       sqrt(numpy.linalg.det(random_mom2[4:6,4:6])))
f.close()
f = open("test_diagnostics_get_emitxy.icc","w")
f.write("BOOST_CHECK_CLOSE(diagnostics.get_emitxy(), %0.15f, tolerance_emit4d);\n" % \
       sqrt(numpy.linalg.det(random_mom2[0:4,0:4])))
f.close()
f = open("test_diagnostics_get_emitxyz.icc","w")
f.write("BOOST_CHECK_CLOSE(diagnostics.get_emitxyz(), %0.15f, tolerance_emit6d);\n" % \
       sqrt(numpy.linalg.det(random_mom2)))
f.close()
