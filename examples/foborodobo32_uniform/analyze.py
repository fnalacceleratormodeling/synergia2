#!/usr/bin/env python3
import sys, os

import numpy as np
import matplotlib.pyplot as plt
import h5py

h5d = h5py.File('diag_b000.h5', 'r')
mean = h5d.get('mean')[()]
std = h5d.get('std')[()]
mins = h5d.get('min')[()]
maxs = h5d.get('max')[()]
numpart = h5d.get('num_particles')[()]
s = h5d.get('s')[()]

mass = h5d.get('mass')
pz = h5d.get('pz')[()]

betagamma = pz[0]/mass
gamma = np.sqrt(betagamma**2 + 1)
beta = betagamma/gamma

print('initial pz: ', pz[0])
print('gamma: ', gamma)
print('beta: ', beta)

print('initial macroparticles: ', numpart[0])
print('turns: ', s.shape[0]-1)

plt.figure()
plt.title('transverse means')
plt.plot(mean[:,0], label='<x>')
plt.plot(mean[:,2], label='<y>')
plt.xlabel('turn')
plt.legend(loc='best')

plt.figure()
plt.title('transverse stds')
plt.plot(std[:,0], label='std(x)')
plt.plot(std[:,2], label='std(y)')
plt.xlabel('turn')
plt.legend(loc='best')

plt.figure()
plt.title('transverse extents')
plt.plot(mins[:,0], label='min x')
plt.plot(mins[:,1], label='min y')
plt.plot(maxs[:,0], label='max x')
plt.plot(maxs[:,1], label='max y')
plt.xlabel('turn')
plt.legend(loc='best')

plt.figure()
plt.title('longitudinal extents')
plt.plot(mins[:,2]*beta, label='min s')
plt.plot(maxs[:,2]*beta, label='max s')
plt.xlabel('turn')
plt.ylabel('s [m]')
plt.legend(loc='best')

turns = list(range(0, s.shape[0]+1, 10))
if turns[-1] != s.shape[0]-1:
    s.append(s.shape[0]-1)

for fnum in [0,1,4,7, 10]:
    h5p = h5py.File('particles_b000_%04d.h5'%fnum, 'r')

    p = h5p.get('particles')[()]
    turn = h5p.get('repetition')[()]

    print('turn ', turn, ' min s: ', p[:,4].min(), ' ,max s: ', p[:,4].max())
    # check for NaN x
    print('x NaN values: ', end='')
    firstnan = True
    for i in range(p.shape[0]):
        if np.isnan(p[i,0]):
            if not firstnan:
                print(', ', end='')
            firstnan = False
            print(i, end='')
    print()
    # check for NaN y
    print('y NaN values: ', end='')
    firstnan = True
    for i in range(p.shape[0]):
        if np.isnan(p[i,2]):
            if not firstnan:
                print(', ', end='')
            firstnan = False
            print(i, end='')
    print()
    # check for NaN cdt
    print('cdt NaN values: ', end='')
    firstnan = True
    for i in range(p.shape[0]):
        if np.isnan(p[i,4]):
            if not firstnan:
                print(', ', end='')
            firstnan = False
            print(i, end='')
    print()
    plt.figure()
    plt.title('longitudinal coord turn %d'%turn)
    plt.hist(p[:,4], 40, label='particle s')
    plt.legend(loc='best')

    plt.figure()
    plt.title('longitudinal dp/p turn %d'%turn)
    plt.hist(p[:,5], 40, label='particle dp/p')
    plt.legend(loc='best')

    h5p.close()

# for turn in turns:
#     h5p = h5py.File('particles_b000_{:04d}.h5'.format(turn), 'r')

#     p = h5p.get('particles')[()]
#     print('turn {:04d} x min: '.format(turn), p[:,0].min())
#     print('turn {:04d} x max: '.format(turn), p[:,0].max())
#     print('turn {:04d} y min: '.format(turn), p[:,2].min())
#     print('turn {:04d} y max: '.format(turn), p[:,2].max())
#     print('turn {:04d} s min: '.format(turn), p[:,4].min()*beta)
#     print('turn {:04d} s max: '.format(turn), p[:,4].max()*beta)
#     del p
#     h5p.close()


h5t = h5py.File('tracks_b000.h5', 'r')
tracks = h5t.get('track_coords')[()]
# tracks file has single particle track data [turn, particle, coord]

plt.figure()
plt.subplot(211)
plt.plot(tracks[:, 0, 0], label='particle 0 x')
plt.legend(loc='best')
plt.subplot(212)
plt.plot(tracks[:, 0, 2], label='particle 0 y')
plt.legend(loc='best')
plt.xlabel('turn')

plt.show()
