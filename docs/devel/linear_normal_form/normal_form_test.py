#!/usr/bin/env python

import numpy as np
from numpy import sin,cos,pi,sqrt
from StringIO import StringIO
import matplotlib.pyplot as plt
import sys
import os

import linear_normal_form

mi_map = np.array([
    [  -6.37425967e-01,   4.73275047e+00,   0.00000000e+00,   0.00000000e+00,   7.13941811e-06,   1.14178418e-01],
    [ -5.69292251e-02,  -1.14612241e+00,   0.00000000e+00,   0.00000000e+00,  -1.43012713e-07,  -2.28715635e-03],
    [ -2.55313450e-16 ,  2.28411626e-15,  -2.12883095e+00,   3.02765363e+01,   4.60831942e-22,   1.42703583e-16],
    [  9.06223330e-18,   2.32731536e-16,  -6.15504874e-02,   4.05638393e-01,   5.82178883e-22,   1.72385663e-18],
    [  8.00209457e-03,   1.20703300e-01,  -9.59798251e-17,  -2.14138188e-15,   9.98154723e-01,  -2.95344987e+01],
    [  5.00359134e-07,   7.54739876e-06,   0.00000000e+00,   0.00000000e+00,   1.24941659e-04,   9.98151777e-01]], 'd')

pathologic_map = np.array([
    [  1.87820723e+00,  -3.01328479e+01,   0.00000000e+00,   0.00000000e+00,  -9.06053010e-02,   6.84168359e+00],
    [  1.85358632e-01,  -2.43052756e+00,   0.00000000e+00,   0.00000000e+00,  -6.80491154e-03,   2.90492489e-01],
    [  0.00000000e+00,   0.00000000e+00,   1.02216594e+00,   3.18724308e-01,   0.00000000e+00,   0.00000000e+00],
    [  0.00000000e+00,   0.00000000e+00,  -1.10126474e-02,   9.74880850e-01,   0.00000000e+00,   0.00000000e+00],
    [ -6.90388840e-01,   6.84265746e+00,   0.00000000e+00,   0.00000000e+00,   8.34990419e-01,   2.42247521e+01],
    [  6.20835190e-03,  -9.11709903e-02,   0.00000000e+00,   0.00000000e+00,  -1.33542083e-02,   8.34556466e-01]], 'd')


remap = linear_normal_form.sensible_order(mi_map)

# compute normal form matrices
(b, binv) = linear_normal_form.normal_form(remap)

#print "B matrix"
#print np.array2string(b, max_line_width=200)


# Generate distributions

# normalized emittance of 18 (pi) mm-mr (95%)
# approximate lattice functions and beam parameters for beam width calculations
betagamma = 9.4
emit = 18.0e-6/(betagamma*6.0*np.pi)
beta_x = 10
beta_y = 60

sigx = np.sqrt(emit*beta_x)
sigy = np.sqrt(emit*beta_y)
sigz = 0.3

moments = np.array([sigx**2,sigy**2, sigz**2],'d')

# generate mean actions for normal form generation
actions = linear_normal_form.get_mean_actions(b, moments)
print "actions: ", actions

# see if any are <= 0
problems = map(lambda(x):x<=0.0, actions).count(True)
if problems != 0:
    raise RuntimeError, "A matching beam distribution cannot be constructed with the given moments."
    

npart = 100000
nturns=100
# generate particles
particles = linear_normal_form.generate_particles(npart, actions, b)

# these particles are generated in "sensible-ordering" (Michelotti)

print "sigx should be: ", sigx, " actual: ", np.std(particles[0,:])
print "sigy should be: ", sigy, " actual: ", np.std(particles[1,:])
print "sigz should be: ", sigz, " actual: ", np.std(particles[2,:])


# map particles around turns keep statistics
meanx =np.zeros((nturns),'d')
meany =np.zeros((nturns),'d')
stdx = np.zeros((nturns),'d')
stdy = np.zeros((nturns),'d')
cxxp = np.zeros((nturns),'d')
cxy = np.zeros((nturns),'d')
cyyp = np.zeros((nturns),'d')
# track action and angle for 10 particles

f1 = plt.figure()
plt.title("turns 0-40 phase space x-x'")
f2 = plt.figure()
plt.title("turns 0-40 phase space y-y''")
f3 = plt.figure()
plt.title("turns 50-90 phase space x-x'")
f4 = plt.figure()
plt.title("turns 50-90 phase space y-y''")

sample_particles = range(20,1000,10)
for turn in range(nturns):
    c = np.cov(particles)
    meanx[turn] = particles[0,:].mean()
    meany[turn] = particles[1,:].mean()
    stdx[turn] = c[0,0]
    stdy[turn] = c[1,1]
    cxxp[turn] = c[0,2]/sqrt(c[0,0]*c[2,2])
    cxy[turn] =  c[0,1]/sqrt(c[0,0]*c[1,1])
    cyyp[turn] = c[1,3]/sqrt(c[1,1]*c[3,3])

    if turn%10 == 0:
        if turn/10 < 5:
            plt.figure(1)
        else:
            plt.figure(3)
        plt.plot(particles[0,:],particles[2,:],'.')
        if turn/10 < 5:
            plt.figure(2)
        else:
            plt.figure(4)
        plt.plot(particles[1,:],particles[3,:],'.')
        
    particles = np.dot(remap, particles)


plt.figure()
plt.plot(stdx, label='std x')
plt.plot(stdy, label='std y')
plt.legend()

plt.figure()
plt.plot(cxxp, label='x-xp corr')
plt.plot(cyyp, label='y-yp corr')
plt.legend()

plt.figure()
plt.plot(cxy, label='x-y corr')
plt.legend()

#plt.figure()
#plt.plot(exes, label='xcoordinate particle 0')
#plt.plot(wyes, label='ycoordinate particle 1')
#plt.legend()
#plt.figure()

for j in []: # range(len(sample_particles)):
    plt.figure()
    #plt.polar(angparticles[0,:,sample_particles[j]],actparticles[0,:,sample_particles[j]])
    plt.polar(angparticles[0,:,j],actparticles[0,:,j],'.')
    plt.polar(angparticles[1,:,j],actparticles[1,:,j],'.')

#plt.figure()
#plt.plot(meanx, label='mean x')
#plt.plot(meany, label='mean y')

plt.show()

