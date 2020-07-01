#!/usr/bin/env python

import os
import sys
import re
import glob
import mpi4py
import mpi4py.MPI as MPI
import tables
import numpy as np
from  tune_suite import *

DEBUG=0

# return type of file (single track or multitrack)
def track_filetype(hfname):
    partid_re = re.compile(basefile + "_([0-9]+).h5")
    partid_mo = partid_re.match(hfname)
    if not partid_mo:
        raise RuntimeError, "filename regular expression match did not include particle ID or rank"
    h5f = tables.openFile(hfname)
    if "coords" in dir(h5f.root):
        # this is a single track file.
        h5f.close()
        return partid_mo.group(1)
    elif "track_coords" in dir(h5f.root):
        # this is a multi track file
        ntracks = h5f.root.track_coords.shape[0]
        h5f.close()
        return -ntracks
    else:
        h5f.close()
        raise RuntimeError, "track_filetype couldn't determine type of file %s"%hfname

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise RuntimeError,"usage: calc_tunes.py tracks-filename"

    os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

    commsize = MPI.COMM_WORLD.Get_size()
    myrank = MPI.COMM_WORLD.Get_rank()

    tracks_file = sys.argv[1]
    if DEBUG and myrank == 0:
        print "Reading tracks from file ", tracks_file

    tunes_file = os.path.splitext(tracks_file)[0]+"_tunes"
    h5f = tables.openFile(tracks_file)

    # collect tunes in a dictionary indexed by particle ID containing
    # a tuple of ([starting turn], [xtunes], [ytunes], [ltunes])

    ntracks = h5f.root.track_coords.shape[0]
    if DEBUG: print "rank: ", myrank, ", ntracks: ", ntracks

    tunelist = {}
    # divide track processing up by processor
    tracks_per_proc = (ntracks+commsize-1)/commsize
    if DEBUG and myrank == 0: print "tracks per proc: ", tracks_per_proc

    my_first_track = myrank*tracks_per_proc
    my_last_track = min( (myrank+1)*tracks_per_proc, ntracks )
    if DEBUG>1: print "proc: ", myrank,", first track: ", my_first_track,", last  track: ", my_last_track

    window_size = 1000
    for do_track in range(my_first_track, my_last_track):
        if DEBUG>1: print "proc: ", myrank,", working on track: ", do_track

        # this is a file of bulk tracks.  This will contain an array
        #[trknum, coords, turn_number] with shape
        # ntracks x 7 x nturns

        coords = h5f.root.track_coords[do_track,0:6,:]
        nturns = coords.shape[1]
        tune_starts = range(0, nturns-windowsize+1, window_size/20)
        xtunes = []
        ytunes = []
        ltunes = []
        for tstart in tune_starts:
            if DEBUG>2:
                print "proc: ", myrank, ", track: ", do_track, ", turn: ", tstart
            #tunes = interp_tunes(coords[:, tstart:tstart+window_size])
            tunes = interp_tunes(coords[:, tstart:tstart+window_size])
            xtunes.append(tunes[0])
            ytunes.append(tunes[1])
            ltunes.append(tunes[2])

        tunelist[do_track] = (tune_starts, xtunes, ytunes, ltunes)
        #tunelist[trknum] = basic_tunes(coords)
        #tunelist[trknum] = cft_tunes(coords)
        if DEBUG>1: print "tunes for particle ", do_track,": ", tunelist[do_track]
                                                                    
    h5f.close()
    # send my tunes to rank 0 for writing out
 
    # rank 0 will collect all the tune data and write it out
    if myrank == 0:
        for r in range(1,commsize):
            rtunes = MPI.COMM_WORLD.recv(source=r)
            if DEBUG: print "Receiving tune data for %d tracks from rank %d"%(len(rtunes),r)
            tunelist.update(rtunes)

        if myrank == 0: print "Saving data for ",len(tunelist), " particles"
        np.save(tunes_file, tunelist)
    else:
        # send my data to rank 0
        MPI.COMM_WORLD.send(tunelist,dest=0)
