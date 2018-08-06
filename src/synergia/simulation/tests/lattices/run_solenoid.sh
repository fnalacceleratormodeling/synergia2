#!/bin/bash

madx madx_track_solenoid.madx
python convert_madx_solenoid_tracks_to_np.py track_solenoidone
