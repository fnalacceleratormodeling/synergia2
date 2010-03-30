#include "propagator.h"

void
Propagator::construct()
{
    Lattice_element_slices all_slices;
    std::cout << "jfa: start triple loop\n";
    for (Steps::const_iterator s_it = stepper.get_steps().begin(); s_it
            != stepper.get_steps().end(); ++s_it) {
        std::cout << "jfa: step\n";

        for (Operators::const_iterator o_it = (*s_it)->get_operators().begin(); o_it
                != (*s_it)->get_operators().end(); ++o_it) {
            std::cout << "jfa: operator " << (*o_it)->get_name() << "\n";

            for(Lattice_element_slices::const_iterator les_it = (*o_it)->get_slices().begin();
                    les_it != (*o_it)->get_slices().end(); ++les_it) {
                std::cout << "jfa: add slice\n";
                all_slices.push_back(*les_it);
            }
        }
    }
    std::cout << "jfa: end triple loop\n";
    chef_lattice.construct_sliced_beamline(all_slices);
}

Propagator::Propagator(Stepper & stepper) :
    stepper(stepper), chef_lattice(stepper.get_lattice())
{
    construct();
}

void
Propagator::propagate(Bunch & bunch, int num_turns, bool diagnostics_per_step,
        bool diagnostics_per_turn)
{

}

Chef_lattice &
Propagator::get_chef_lattice()
{
    return chef_lattice;
}

Propagator::~Propagator()
{

}
