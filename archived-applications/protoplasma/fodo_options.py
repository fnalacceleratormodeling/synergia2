#!/usr/bin/env python

from math import pi
import synergia_workflow

opts = synergia_workflow.Options("protoplasma")

opts.add("seed", 4, "Pseudorandom number generator seed", int)
opts.add("radius", 0.05, "aperture radius [m]", float)
opts.add("real_particles", 1.0e11, "Number of real particles", float)
opts.add("macro_particles", 100000, "Number of macro particles", int)
opts.add("verbose", False, "Verbose propagation", bool)

# diagnostics saving options
opts.add("step_tracks", 0, "Number of particles to track each step", int)
opts.add("step_full2", True, "Whether to do full2 diagnostics each step", bool)
opts.add("step_particles", False, "Whether to save all particles each step", bool)
opts.add("turn_tracks", 0, "Number of particles to track each turn", int)
opts.add("turn_full2", True, "Whether to do full2 diagnostics each turn", bool)
opts.add("turn_particles", False, "Whether to save all particles each turn", bool)

#opts.add("num_steps", 728, "Number of steps per turn", int)
opts.add("num_steps", 120, "Number of steps per turn", int)
opts.add("num_turns", 1, "Number of turns", int)
opts.add("map_order", 4, "Map order", int)
opts.add("partpercell", 10, "macro particles per grid cell", int)
opts.add("lattice_load", False, "load lattice setting", bool)

opts.add("steps", 240, "Number of steps per turn", int)
opts.add("turns", 1, "Number of turns", int)
opts.add("emit", 1e-6, "Horizontal and vertical emittance [m rad]", float)
opts.add("stdz", 0.01, "RMS longitudinal length [m]", float)
opts.add("dpop", 1e-4, "delta p / p", float)

opts.add("deltae", 1.0, "energy offset", float)

opts.add("sigx", 20e-6, "Horizontal beam size [m]", float)
opts.add("sigy", 20e-6, "Vertical beam size [m]", float)
opts.add("bunchlen", 4.6, "RMS bunch length (z width) [nanoseconds]", float)

opts.add("x_offset", 0.0, "Bunch offset in x", float)
opts.add("y_offset", 0.0, "Bunch offset in y", float)
opts.add("z_offset", 0.0, "Bunch offset in z", float)
opts.add("xp_offset", 0.0, "Bunch offset in xp", float)
opts.add("yp_offset", 0.0, "Bunch offset in yp", float)
opts.add("zp_offset", 0.0, "Bunch offset in zp", float)

#opts.add("spacecharge", True, "Use hockney 3d open space charge", bool)
opts.add("spacecharge", "no_op", "Use hockney 3d open space charge", str) 

opts.add("gridx", 32, "size of transverse grid for solver", int)
opts.add("gridy", 32, "size of transverse grid for solver", int)
opts.add("gridz", 128, "size of longitudinal grid for solver", int)
opts.add("save_tracks", 0, "Number of particles to track each turn", int)

job_mgr = synergia_workflow.Job_manager("fodo.py", opts, ["tevatron.lat"])
