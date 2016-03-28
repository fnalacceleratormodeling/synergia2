#!/bin/bash

seqs="mpole_k1 mpole_k1s mpole_k1_tilt mpole_k2 mpole_k2s mpole_k2_tilt mpole_k3 mpole_k3s mpole_k3_test"

for s in $seqs
do
    cp $s.seq mpole1.seq
    madx madx_track_multipole.madx
    python convert_madx_tracks_to_np.py madx_track_multipole.txtone
    mv madx_track_multipole.npy $s.npy
done
