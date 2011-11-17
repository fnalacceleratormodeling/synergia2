#!/usr/bin/env python

from math import pi
import synergia_workflow

opts = synergia_workflow.Options("mu2e")

opts.add("seed", 4, "Pseudorandom number generator seed", int)
opts.add("radius", 0.05, "aperture radius [m]", float)
opts.add("real_particles", 3.0e12, "Number of real particles", float)
opts.add("macro_particles", 0, "Number of macro particles", int)
opts.add("rampturns", 100, "sextupole ramping turns", int)
opts.add("extraction_fraction", 7.0/8.0, "extraction fraction", float)
opts.add("verbose", False, "Verbose propagation", bool)

# diagnostics saving options
opts.add("step_tracks", 0, "Number of particles to track each step", int)
opts.add("step_full2", True, "Whether to do full2 diagnostics each step", bool)
opts.add("step_particles", False, "Whether to save all particles each step", bool)
opts.add("turn_tracks", 0, "Number of particles to track each turn", int)
opts.add("turn_full2", True, "Whether to do full2 diagnostics each turn", bool)
opts.add("turn_particles", False, "Whether to save all particles each turn", bool)
opts.add("twiss", False, "Whether to save twiss parameters each turn", bool)
opts.add("separatrix", False, "Whether to save SFP points of separatrix lines each turn", bool)

#opts.add("num_steps", 728, "Number of steps per turn", int)
opts.add("num_steps", 2, "Number of steps per turn", int)
opts.add("num_turns", 1, "Number of turns", int)
opts.add("map_order", 2, "Map order", int)
opts.add("partpercell", 10, "macro particles per grid cell", int)

# normalized geometric emittance of 19.09e-6 mm-mr at beta*gamma=9.4855
# gives a beam spot with a half width of 14 mm at the ipm location with beta_x
# is 51.67.
opts.add("norm_emit",20e-6, "Horizontal and vertical emittance [m rad]", float)
opts.add("bunchlen", 40.0, "RMS bunch length (z width) [nanoseconds]", float)
opts.add("stdz", 0.3, "RMS longitudinal length [m]", float)
opts.add("tuneh", 9.65, "horizontal fractional tune", float)
opts.add("tunev", 9.75, "vertical fractional tune", float)
opts.add("resonant_tune", 29.0/3.0, "resonant tune", float)
opts.add("phase_g", 300, "phase of g in deg", float)
opts.add("flip", False, "Whether to change polarity of harmonic circuits", bool)

# RF voltage is 32 kV divided among 3 cavities
# rf cavity voltage, is 32.0 kV total distributed over 3 cavities. 
# MAD8 expects cavities voltages in  units of MV.
opts.add("rf_voltage", 0.0, "RF cavity voltage in MV", float)

opts.add("x_offset", 0.0, "Bunch offset in x", float)
opts.add("y_offset", 0.0, "Bunch offset in y", float)
opts.add("z_offset", 0.0, "Bunch offset in z", float)

#opts.add("spacecharge", True, "Use hockney 3d open space charge", bool)
opts.add("spacecharge", "no_op", "Use hockney 3d open space charge", str) 

opts.add("gridx", 32, "size of transverse grid for solver", int)
opts.add("gridy", 32, "size of transverse grid for solver", int)
opts.add("gridz", 128, "size of longitudinal grid for solver", int)
opts.add("save_tracks", 0, "Number of particles to track each turn", int)

job_mgr = synergia_workflow.Job_manager("mu2e.py", opts, ["Debunch_modified_rf.lat"])
