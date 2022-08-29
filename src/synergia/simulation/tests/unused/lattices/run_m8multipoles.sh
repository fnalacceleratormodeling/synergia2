#!/bin/bash

seqs="mpole_k1 mpole_k1_tilt mpole_k2 mpole_k2_tilt mpole_k3 mpole_k3_tilt mpole_k4 mpole_k4_tilt mpole_k5 mpole_k5_tilt mpole_k6 mpole_k6_tilt"

for s in $seqs
do
    cp $s.lat mpole1.lat
    mad8 <mad8_track_multipole.mad
    python convert_madx_tracks_to_np.py mad8_track_multipole.txt
    mv mad8_track_multipole.npy m8$s.npy
done
