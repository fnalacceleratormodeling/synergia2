#!/bin/bash

seqs="m_base_quad m_skew_quad m_tilt_quad"

for s in $seqs
do
    echo "use, $s;" >use_misc.mad
    mad8 < mad8_track_misc_elements.mad
    python convert_madx_tracks_to_np.py mad8_track_misc_element.txt
    mv mad8_track_misc_element.npy m8$s.npy
done
