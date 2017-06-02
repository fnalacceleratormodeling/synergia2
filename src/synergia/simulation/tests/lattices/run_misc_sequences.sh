#!/bin/bash

seqs="m_base_quad m_skew_quad m_tilt_quad"

for s in $seqs
do
    echo "use, sequence=$s;" >use_misc.madx
    madx madx_track_misc_element.madx
    python convert_madx_tracks_to_np.py madx_track_multipole.txtone
    mv madx_track_multipole.npy $s.npy
done

s=m_nllens
echo "use, sequence=m_nllens;" >use_misc.madx
madx madx_thintrack_misc_element.madx
python convert_madx_tracks_to_np.py madx_track_multipole.txtone
mv madx_track_multipole.npy $s.npy
