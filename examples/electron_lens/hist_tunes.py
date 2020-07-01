#!/usr/bin/env python
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise RuntimeError,"usage: hist_tunes.py <tunes-filename>"

    tunedict = np.load(sys.argv[1])[()]

    # the format of this stuff
    # tunedict is dictionary indexed by particle id
    # each entry is a tuple of 4 lists
    # tunedict[pid][0]: list of starting turn numbers
    # tunedict[pid][1]: x tunes starting at corresponding turn number
    # tunedict[pid][2]: same for y
    # tunedict[pid][3]: same for longitudinal

    # getting the first tunes for all the coordinates
    xtunes = np.array([t[1][0]*12 for t in tunedict.values()])
    ytunes = np.array([t[2][0]*12 for t in tunedict.values()])
    ltunes = np.array([t[3][0]*12 for t in tunedict.values()])

    # sometimes you want to set the range explicitly, but if you're
    # not careful you can be burned.
    #plt.hist(xtunes,40,range=[0.40,0.45])
    plt.title("tune histograms")
    plt.subplot(211)
    ax1 = plt.gca()
    plt.hist(xtunes,50, range=[3.5, 4.0])
    #plt.hist(xtunes, 50, range=[0.32,0.46])
    plt.xlabel("x tune")

    plt.subplot(212)
    #plt.hist(ytunes,40,range=[0.40,0.45])
    plt.hist(ytunes,50, range=[3.5, 4.0])
    plt.xlabel("y tune")

    plt.savefig("hist_tunes.png")

    plt.show()

