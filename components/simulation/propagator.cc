#include "propagator.h"

//void
//Propagator::construct()
//{
//    Lattice_element_slices all_slices;
//    for (Steps::const_iterator s_it = stepper.get_steps().begin(); s_it
//            != stepper.get_steps().end(); ++s_it) {
//        for (Operators::const_iterator o_it = (*s_it)->get_operators().begin(); o_it
//                != (*s_it)->get_operators().end(); ++o_it) {
//            all_slices.splice(all_slices.end(), (*o_it)->get_slices());
//        }
//    }
//    chef_lattice.construct_sliced_beamline(all_slices);
//}

Propagator::Propagator(Stepper_sptr const& stepper_sptr) :
    stepper_sptr(stepper_sptr)
{
}

void
Propagator::propagate(Bunch & bunch, int num_turns, bool diagnostics_per_step,
        bool diagnostics_per_turn)
{
    for (int turn = 0; turn < num_turns; ++turn) {
        for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
                != stepper_sptr->get_steps().end(); ++it) {
            std::cout << "jfa: step\n";
            (*it)->apply(bunch);
            if (diagnostics_per_step) {
                std::cout << "jfa: per-step diagnostics\n";
            }
        }
        if (diagnostics_per_turn) {
            std::cout << "jfa: per-turn diagnostics\n";
        }
    }
}

Propagator::~Propagator()
{

}
