#!/usr/bin/env python
import synergia
from cold_fodo_options import opts
import sys

lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)
bunch = synergia.optics.generate_matched_bunch_transverse(
              lattice_simulator, opts.emit, opts.emit, opts.stdz, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)
# make the beam cold in the transverse dimensions
bunch.get_local_particles()[:,1] = 1.0e-10
bunch.get_local_particles()[:,3] = 1.0e-10
bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
for element in lattice.get_elements():
    if opts.aperture == "circular":
        element.set_string_attribute("aperture_type","circular")
        element.set_double_attribute("circular_aperture_radius", 0.005)
    elif opts.aperture == "elliptical":
        element.set_string_attribute("aperture_type","elliptical")
        element.set_double_attribute("elliptical_aperture_horizontal_radius", 0.005)
        element.set_double_attribute("elliptical_aperture_vertical_radius", 0.002)
    elif opts.aperture == "rectangular":
        element.set_string_attribute("aperture_type","rectangular")
        element.set_double_attribute("rectangular_aperture_width", 2*0.005)
        element.set_double_attribute("rectangular_aperture_height", 2*0.002)
    elif opts.aperture == "polygon":
        element.set_string_attribute("aperture_type","polygon")
        element.set_double_attribute("the_number_of_vertices", 4)

        element.set_double_attribute("pax1", 0.005)
        element.set_double_attribute("pay1", 0.0)
        element.set_double_attribute("pax2", 0.0)
        element.set_double_attribute("pay2", 0.005)
        element.set_double_attribute("pax3", -0.005)
        element.set_double_attribute("pay3", 0.0)
        element.set_double_attribute("pax4", 0.0)
        element.set_double_attribute("pay4", -0.005)
    else:
        print "unknown aperture type '%s'" % opts.aperture
        sys.exit(1)

# turn off magnets
for element in lattice.get_elements():
    if element.has_double_attribute("k1"):
        element.set_double_attribute("k1", 0.0)
lattice_simulator.update()

stepper = synergia.simulation.Independent_stepper_elements(
                            lattice_simulator, opts.steps)
diagnostics_step = synergia.bunch.Diagnostics_full2("full2.h5")
diagnostics_turn = synergia.bunch.Diagnostics_particles("particles.h5")
bunch_simulator.add_per_step(diagnostics_step)
bunch_simulator.add_per_turn(diagnostics_turn)
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch_simulator, opts.turns, opts.turns,
                     opts.verbose)
