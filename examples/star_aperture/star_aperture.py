#!/usr/bin/env python
import synergia
from star_aperture_options import opts
from vertex import vertex
import sys
import numpy as np

lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)
bunch = synergia.optics.generate_matched_bunch_transverse(
              lattice_simulator, opts.emit, opts.emit, opts.stdz, opts.dpop,
              opts.real_particles, opts.macro_particles,
              seed=opts.seed)
# make the beam cold in the transverse dimensions
bunch.get_local_particles()[:,1] = 1.0e-10
bunch.get_local_particles()[:,3] = 1.0e-10

circular_aperture = 3.2848 / 2 * 0.0254
wire_x = 0.016
wire_width = 0.0001
gap = 0.014

x = np.zeros(20,'d')
y = np.zeros(20,'d')
(x, y) = vertex()
u = 2.811 * 0.0254
v = 0.0

for element in lattice.get_elements():
    if element.get_type() == "drift":
        element.set_string_attribute("aperture_type","circular")
        element.set_double_attribute("circular_aperture_radius", circular_aperture)
    elif element.get_type() == "quadrupole":
        element.set_string_attribute("aperture_type","polygon")
        element.set_double_attribute("the_number_of_vertices", 72)

        element.set_double_attribute("pax1", u)
        element.set_double_attribute("pay1", v)

        for i in range(0, 17):
            pax = "pax" + str(i + 17 * 0 + 2)
            pay = "pay" + str(i + 17 * 0 + 2)
            element.set_double_attribute(pax, x[i])
            element.set_double_attribute(pay, y[i])

        element.set_double_attribute("pax19", v)
        element.set_double_attribute("pay19", u)

        for i in range(0, 17):
            pax = "pax" + str(i + 17 * 1 + 3)
            pay = "pay" + str(i + 17 * 1 + 3)
            element.set_double_attribute(pax, -x[16-i])
            element.set_double_attribute(pay, y[16-i])

        element.set_double_attribute("pax37", -u)
        element.set_double_attribute("pay37", v)

        for i in range(0, 17):
            pax = "pax" + str(i + 17 * 2 + 4)
            pay = "pay" + str(i + 17 * 2 + 4)
            element.set_double_attribute(pax, -x[i])
            element.set_double_attribute(pay, -y[i])

        element.set_double_attribute("pax55", v)
        element.set_double_attribute("pay55", -u)

        for i in range(0, 17):
            pax = "pax" + str(i + 17 * 3 + 5)
            pay = "pay" + str(i + 17 * 3 + 5)
            element.set_double_attribute(pax, x[16-i])
            element.set_double_attribute(pay, -y[16-i])
    elif element.get_name() == "e_septum":
        element.set_string_attribute("aperture_type","wire_elliptical")
        element.set_double_attribute("wire_elliptical_aperture_horizontal_radius", circular_aperture)
        element.set_double_attribute("wire_elliptical_aperture_vertical_radius", circular_aperture)
        element.set_double_attribute("wire_elliptical_aperture_wire_x", wire_x)
        element.set_double_attribute("wire_elliptical_aperture_wire_width", wire_width)
        element.set_double_attribute("wire_elliptical_aperture_gap", gap)
    #else:
    #    print "unknown element '%s'" % element.get_name()
    #    sys.exit(1)

    lattice.print_()

# turn off magnets
for element in lattice.get_elements():
    if element.has_double_attribute("k1"):
        element.set_double_attribute("k1", 0.0)
lattice_simulator.update()
bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

stepper = synergia.simulation.Independent_stepper_elements(
                            lattice_simulator, opts.steps)
diagnostics_step = synergia.bunch.Diagnostics_full2("full2.h5")
diagnostics_turn = synergia.bunch.Diagnostics_particles("particles.h5")
bunch_simulator.add_per_step(diagnostics_step)
bunch_simulator.add_per_turn(diagnostics_turn)
propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch_simulator, opts.turns, opts.turns, opts.verbose)
