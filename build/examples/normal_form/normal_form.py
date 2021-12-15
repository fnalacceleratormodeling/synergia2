

import numpy as np
import synergia as syn

t = syn.foundation.Trigon_o2(2.0)
print(f't.power = {t.power()}, t.dim = {t.dim()}, t.term_count = {t.count()}')
print(t.value())
t.set(3.0)
print(t.value())

print(t.get_term(0, 0))
t.set_term(0, 0, 4.0)
print(t.get_term(0, 0))

print(syn.foundation.Trigon_index_to_canonical_o4(12))
print(syn.foundation.Trigon_canonical_to_index([0, 0, 3, 2]))

print(t.get_term([0, 2]))

m = syn.foundation.TMapping_o2()
m.component(0).print_coeffs()


def run():

    # logger
    screen = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.DEBUG)

    simlog = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.INFO_TURN)

    # lattice
    reader = syn.lattice.MadX_reader()
    lattice = reader.get_lattice("fodo", "channel.madx")
    lattice.set_all_string_attribute("extractor_type", "libff")
    syn.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    # normal form
    nf = syn.simulation.Lattice_simulator.calculate_normal_form_o2(lattice)

    #nf.get_f()[0].component(0).print_coeffs()
    #nf.get_g()[0].component(0).print_coeffs()

    ref = lattice.get_reference_particle()
    energy = ref.get_total_energy()
    momentum = ref.get_momentum()
    gamma = ref.get_gamma()
    beta = ref.get_beta()

    macro_particles = 80
    real_particles = 1.0e9

    stepper = syn.simulation.Independent_stepper_elements()
    propagator = syn.simulation.Propagator(lattice, stepper)

    sim = syn.simulation.Bunch_simulator.create_single_bunch_simulator( 
            ref, macro_particles, real_particles)


    bunch = sim.get_bunch()
    hparts = bunch.get_particles_numpy()

    offset = 0.001

    for i in range(macro_particles//2):
        hparts[i, 0:6] = 0.0
        hparts[i, 0] = offset * i

    for i in range(macro_particles//2, macro_particles):
        hparts[i, 0:6] = 0.0
        hparts[i, 2] = offset*(i-macro_particles/2)

    bunch.checkin_particles()

    bunch.print_particle(0, screen)
    bunch.print_particle(1, screen)
    bunch.print_particle(2, screen)
    bunch.print_particle(3, screen)

    bunch.print_particle(macro_particles//2 + 0, screen)
    bunch.print_particle(macro_particles//2 + 1, screen)
    bunch.print_particle(macro_particles//2 + 2, screen)
    bunch.print_particle(macro_particles//2 + 3, screen)

    fh = open("nf.dat", "w")

    def action(sim, lattice, turn):
        bunch = sim.get_bunch()
        bunch.checkout_particles()
        hparts = bunch.get_particles_numpy()

        for p in range(bunch.size()):
            hp = np.zeros(6, dtype='d')
            hp[:] = hparts[p, 0:6]
            tm = nf.convert_xyz_to_normal(hp)

            print('{0} {1} {2} {3} {4} {5} '.format(
                tm[0].real, tm[0].imag,
                tm[1].real, tm[1].imag,
                tm[2].real, tm[2].imag), file=fh, end="")

        print(file=fh)
        fh.flush()

    # register turn end action for writing particles out
    sim.reg_prop_action_turn_end(action)

    # propagate
    propagator.propagate(sim, simlog, 10);

def main():
    run()

main()



