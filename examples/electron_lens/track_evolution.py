#!/usr/bin/env python
import sys,os
import numpy as np
import tables
import matplotlib.pyplot as plt

import tune_suite

os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

def calculate_actions(coords):
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
def calculate_emittances(p):
    beta_x = 17.3259852015
    alpha_x = 1.85063532729
    beta_y = 17.2246414528
    alpha_y = -1.90226103118

    xactions = (p[:,0]**2 + (alpha_x*p[:,0]+beta_x*p[:,1])**2)/(2.0*beta_x)
    yactions = (p[:,2]**2 + (alpha_y*p[:,2]+beta_y*p[:,3])**2)/(2.0*beta_y)

    npart = p.shape[0]

    xactions.sort()
    yactions.sort()
    i95 = int(npart * 0.95 + 1)
    i99 = int(npart * 0.99 + 1)
    i999 = int(npart * 0.999 + 1)
    xemitrms = xactions.mean()
    xemit95 = xactions[i95]
    xemit99 = xactions[i99]
    xemit999 = xactions[i999]
    yemitrms = yactions.mean()
    yemit95 = yactions[i95]
    yemit99 = yactions[i99]
    yemit999 = yactions[i999]
    return (xemitrms, xemit95, xemit99, xemit999, yemitrms, yemit95, yemit99, yemit999)

if __name__ == "__main__":
    os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

    hfile_nosc = "/data/egstern/elens/elens_nosc.00"
    hfile_3p33e10 = "/data/egstern/elens/elens_3.33e10.00"
    hfile_1pe11 = "/data/egstern/elens/elens_1.0e11.00"

    h5_nosc = tables.openFile(hfile_nosc+"/bulk_track_0.h5")
    h5_3p33e10 = tables.openFile(hfile_3p33e10+"/bulk_track_0.h5")
    h5_1pe11 = tables.openFile(hfile_1pe11+"/bulk_track_0.h5")

    # particle indices to track
    lowxidx = 2782
    medxidx = 2085
    highxidx = 2884

    lowyidx = 1349
    medyidx = 2761
    highyidx = 7142

    # the actions are in order [[lowx, medx, highx],[lowy, medy, highy]]
    actions = np.zeros((2,3), dtype='d')
    actions[0,0] = np.array(calculate_actions(h5_nosc.root.track_coords[lowxidx,:,0]))[0]
    actions[0,1] = np.array(calculate_actions(h5_nosc.root.track_coords[medxidx,:,0]))[0]
    actions[0,2] = np.array(calculate_actions(h5_nosc.root.track_coords[highxidx,:,0]))[0]

    actions[1,0] = np.array(calculate_actions(h5_nosc.root.track_coords[lowyidx,:,0]))[1]
    actions[1,1] = np.array(calculate_actions(h5_nosc.root.track_coords[medyidx,:,0]))[1]
    actions[1,2] = np.array(calculate_actions(h5_nosc.root.track_coords[highyidx,:,0]))[1]

    print "x particle actions: ", actions[0, :]
    print "y particle actions: ", actions[1, :]

    plt.figure()
    plt.title("low x emittance particles")
    plt.plot(h5_1pe11.root.track_coords[lowxidx,0, :], h5_1pe11.root.track_coords[lowxidx, 1, :], '.', label="1.0e11")
    plt.plot(h5_3p33e10.root.track_coords[lowxidx,0, :], h5_3p33e10.root.track_coords[lowxidx, 1, :], '.', label="3.33e10")
    plt.plot(h5_nosc.root.track_coords[lowxidx,0,:], h5_nosc.root.track_coords[lowxidx, 1, :], '.', label="no SC")
    plt.legend(loc='upper right')
    plt.xlabel("X")
    plt.ylabel("XP")

    spchgtrks = [h5_nosc, h5_3p33e10, h5_1pe11]
    xactionidx = [lowxidx, medxidx, highxidx]

    # array of tunes [actions, spcchg]
    spctune = np.zeros((3,3), dtype='d')

    for iaction in range(3):
        for spc in range(3):
            tunes = tune_suite.interp_tunes(spchgtrks[spc].root.track_coords[xactionidx[iaction], :, :])

            #print "iaction: ", iaction, ", spc: ", spc
            #print "spctune: ", spctune.shape
            #print "tunes: ", tunes
            #print "tunes[0]: ", tunes[0]
            #print "tunes[0]*12: ", tunes[0]*12
            spctune[iaction, spc] = tunes[0]*12

    plt.figure()
    plt.title("x tunes vs. action")
    plt.plot(actions[0,:], spctune[:, 0], label='no space charge')
    plt.plot(actions[0,:], spctune[:, 1], label='medium space charge')
    plt.plot(actions[0,:], spctune[:,2], label='large space charge')
    plt.legend(loc='upper right')
    plt.xlabel('X action')
    plt.ylabel('x tune')
    
    bchg = np.array([0.0, 3.33e10, 1.0e11])
    plt.figure()
    plt.title("x tunes vs. bunch charge")
    plt.plot(bchg, spctune[0, :], label='small action')
    plt.plot(bchg, spctune[1, :], label='medium action')
    plt.plot(bchg, spctune[2, :], label='large action')
    plt.legend(loc='upper right')
    plt.xlabel('bunch charge')
    plt.ylabel('x tune')

    print "lowx tunes: ", spctune[0, :]
    print "medx tunes: ", spctune[1, :]
    print "highx tunes: ", spctune[2, :]

    plt.figure()
    plt.title("lowx x emittance particles")
    plt.plot(h5_1pe11.root.track_coords[lowxidx,0, :], h5_1pe11.root.track_coords[lowxidx, 1, :], '.', label="1.0e11")
    plt.plot(h5_3p33e10.root.track_coords[lowxidx,0, :], h5_3p33e10.root.track_coords[lowxidx, 1, :], '.', label="3.33e10")
    plt.plot(h5_nosc.root.track_coords[lowxidx,0,:], h5_nosc.root.track_coords[lowxidx, 1, :], '.', label="no SC")
    plt.legend(loc='upper right')
    plt.xlabel("X")
    plt.ylabel("XP")

    plt.figure()
    plt.title("med x emittance particles")
    plt.plot(h5_1pe11.root.track_coords[medxidx,0, :], h5_1pe11.root.track_coords[medxidx, 1, :], '.', label="1.0e11")
    plt.plot(h5_3p33e10.root.track_coords[medxidx,0, :], h5_3p33e10.root.track_coords[medxidx, 1, :], '.', label="3.33e10")
    plt.plot(h5_nosc.root.track_coords[medxidx,0,:], h5_nosc.root.track_coords[medxidx, 1, :], '.', label="no SC")
    plt.legend(loc='upper right')
    plt.xlabel("X")
    plt.ylabel("XP")

    plt.figure()
    plt.title("high x emittance particles")
    plt.plot(h5_1pe11.root.track_coords[highxidx,0, :], h5_1pe11.root.track_coords[highxidx, 1, :], '.', label="1.0e11")
    plt.plot(h5_3p33e10.root.track_coords[highxidx,0, :], h5_3p33e10.root.track_coords[highxidx, 1, :], '.', label="3.33e10")
    plt.plot(h5_nosc.root.track_coords[highxidx,0,:], h5_nosc.root.track_coords[highxidx, 1, :], '.', label="no SC")
    plt.legend(loc='upper right')
    plt.xlabel("X")
    plt.ylabel("XP")

    plt.show()
