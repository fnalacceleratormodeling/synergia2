#!/usr/bin/env python
import sys,os
import numpy as np
import tables
import matplotlib.pyplot as plt

import tune_suite

os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"
low_emitx = 4.5e-8
med_emitx = 4.25e-6
high_emitx = 8.0e-6

low_emity = 1.1e-7
med_emity = 3.2e-6
high_emity = 6.8e-6

def calculate_action(coords):
    beta_x = 17.3259852015
    alpha_x = 1.85063532729
    beta_y = 17.2246414528
    alpha_y = -1.90226103118

    x = coords[0]
    xp = x*alpha_x + coords[1]*beta_x
    xaction = (x**2 + xp**2)/(2.0*beta_x)
    y = coords[2]
    yp = y*alpha_y + coords[3]*beta_y
    yaction = (y**2 + yp**2)/(2.0*beta_y)
    return (xaction, yaction)


# p is the particle array [npart, coord]
def calculate_actions(p):
    # bpm_beta_x = 4.70792560035
    # bpm_alpha_x = -0.449832076233
    # bpm_beta_y = 46.1594967296
    # bpm_alpha_y = 3.37300523129
    beta_x = 17.3259852015
    alpha_x = 1.85063532729
    beta_y = 17.2246414528
    alpha_y = -1.90226103118

    xactions = (p[:,0]**2 + (alpha_x*p[:,0]+beta_x*p[:,1])**2)/(2.0*beta_x)
    yactions = (p[:,2]**2 + (alpha_y*p[:,2]+beta_y*p[:,3])**2)/(2.0*beta_y)

    return (xactions, yactions)

if __name__ == "__main__":                
    trkfile = "bulk_track_0.h5"
    h5 = tables.openFile(trkfile)
    npart = h5.root.track_coords.shape[0]
    nsteps = h5.root.track_coords.shape[2]

    # select particles from first entry in track file
    particles = h5.root.track_coords[:, :, 0]

    xactions, yactions = calculate_actions(particles)

    minaction = xactions.min()
    minidx = xactions.argmin()
    print "minimum xaction: ", minaction, " at index: ", minidx

    tunes = tune_suite.interp_tunes(h5.root.track_coords[minidx, :, :])
    print "tunes: ", tunes[0]*12, tunes[1]*12

    print
    print

    lowx_list = []
    medx_list = []
    highx_list = []
    lowy_list = []
    medy_list = []
    highy_list = []

    for i, xa in enumerate(xactions):

        if xa > low_emitx * 0.99 and xa < low_emitx * 1.01:
            lowx_list.append(i)

        if xa > med_emitx * 0.99 and xa < med_emitx * 1.01:
            medx_list.append(i)

        if xa > high_emitx * 0.99 and xa < high_emitx * 1.01:
            highx_list.append(i)

    for i, ya in enumerate(yactions):

        if ya > low_emity * 0.99 and ya < low_emity * 1.01:
            lowy_list.append(i)

        if ya > med_emity * 0.99 and ya < med_emity * 1.01:
            medy_list.append(i)

        if ya > high_emity * 0.99 and ya < high_emity * 1.01:
            highy_list.append(i)

    print "low_emitx: ", len(lowx_list), " particles"
    print [(k, particles[k, 6]) for k in lowx_list]
    print "med_emitx: ", len(medx_list), " particles"
    print [(k, particles[k, 6]) for k in medx_list]
    print "high_emitx: ", len(highx_list), " particles"
    print [(k, particles[k, 6]) for k in highx_list]
    print "low_emity: ", len(lowy_list), " particles"
    print [(k, particles[k, 6]) for k in lowy_list]
    print "med_emity: ", len(medy_list), " particles"
    print [(k, particles[k, 6]) for k in medy_list]
    print "high_emity: ", len(highy_list), " particles"
    print [(k, particles[k, 6]) for k in highy_list]
