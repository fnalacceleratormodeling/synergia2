#!/bin/bash

madx madx_track_nllens.madx
python convert_madx_tracks_to_np.py track_nllensone
