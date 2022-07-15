#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("aperture_demo")
opts.add("macroparticles", 20000, "Number of macro particles", int)
opts.add("real_particles", 1.2e12, "Number of real particles", float)

opts.add("aperture", "elliptical", "Aperture type (either 'circular', 'elliptical', 'rectangular', or 'polygon'")

opts.add('hoffset', 0.0025, 'horizontal aperture offset [m]')
opts.add('voffset', 0.001, 'vertical aperture offset [m]')
opts.add('circular_aperture_radius', 0.005, 'circular aperture radius')
opts.add('elliptical_horizontal_radius', 0.005, 'elliptical horizontal radius')
opts.add('elliptical_vertical_radius', 0.002, 'elliptical vertical radius')
opts.add('rectangular_aperture_width', 2*0.005, 'rectangular aperture width')
opts.add('rectangular_aperture_height', 2*0.002, 'rectangular aperture height')
opts.add('polygon_aperture_arm', 0.005, 'length arm of star aperture')

opts.add("seed", 0, "Pseudorandom number generator seed", int)
opts.add("verbose", True, "Verbose propagation", bool)

# no lattice file needed
job_mgr = synergia_workflow.Job_manager("aperture_demo.py", opts, [])

