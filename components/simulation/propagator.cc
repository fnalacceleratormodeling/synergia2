#include "propagator.h"

void
Propagator::construct()
{
    Lattice_element_slices all_slices;
    for (Steps::const_iterator s_it = stepper.get_steps().begin(); s_it
            != stepper.get_steps().end(); ++s_it) {
        for (Operators::const_iterator o_it = (*s_it)->get_operators().begin(); o_it
                != (*s_it)->get_operators().end(); ++o_it) {
            for (Lattice_element_slices::const_iterator les_it =
                    (*o_it)->get_slices().begin(); les_it
                    != (*o_it)->get_slices().end(); ++les_it) {
                all_slices.push_back(*les_it);
            }
        }
    }
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
    for (int turn = 0; turn < num_turns; ++turn) {
        for (Steps::const_iterator step_it = stepper.get_steps().begin(); step_it
                != stepper.get_steps().end(); ++step_it) {
            std::cout << "jfa: step\n";
            for (Operators::const_iterator op_it =
                    (*step_it)->get_operators().begin(); op_it
                    != (*step_it)->get_operators().end(); ++op_it) {
                std::cout << "jfa: apply operator " << (*op_it)->get_name()
                        << std::endl;
                (*op_it)->apply(bunch, chef_lattice);
            }
            if (diagnostics_per_step) {
                std::cout << "jfa: per-step diagnostics\n";
            }
        }
        if (diagnostics_per_turn) {
            std::cout << "jfa: per-turn diagnostics\n";
        }
    }
}

Chef_lattice &
Propagator::get_chef_lattice()
{
    return chef_lattice;
}

Propagator::~Propagator()
{

}
