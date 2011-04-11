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

Propagator::Propagator(Stepper_sptr stepper_sptr) :
    stepper_sptr(stepper_sptr)
{
}

// Object_to_sptr_hack is a local struct
struct Object_to_sptr_hack
{
    void
    operator()(void const *) const
    {
    }
};

void
Propagator::propagate(Bunch & bunch, int num_turns,
        Diagnostics & per_step_diagnostics,
        Diagnostics & per_turn_diagnostics, bool verbose)
{
    Multi_diagnostics multi_per_step_diagnostics;
    multi_per_step_diagnostics.append(Diagnostics_sptr(
            &per_step_diagnostics, Object_to_sptr_hack()));
    Multi_diagnostics multi_per_turn_diagnostics;
    multi_per_turn_diagnostics.append(Diagnostics_sptr(
            &per_turn_diagnostics, Object_to_sptr_hack()));
    propagate(bunch, num_turns, multi_per_step_diagnostics,
            multi_per_turn_diagnostics, verbose);
}

void
Propagator::propagate(Bunch & bunch, int num_turns,
        Multi_diagnostics & per_step_diagnostics,
        Multi_diagnostics & per_turn_diagnostics, bool verbose)
{
    for (int turn = 0; turn < num_turns; ++turn) {
        if (verbose) {
            std::cout << "Propagator: turn " << turn + 1 << "/" << num_turns
                    << std::endl;
        }
        bunch.get_reference_particle().start_repetition();
        for (Multi_diagnostics::iterator dit =
                per_turn_diagnostics.begin(); dit != per_turn_diagnostics.end(); ++dit) {
            std::cout << "jfa&egs: propagator is broken\n";
//            (*dit)->update_and_write(bunch);
        }
        int step_count = 0;
        int num_steps = stepper_sptr->get_steps().size();
        for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
                != stepper_sptr->get_steps().end(); ++it) {
            for (Multi_diagnostics::iterator dit =
                    per_step_diagnostics.begin(); dit
                    != per_step_diagnostics.end(); ++dit) {
                std::cout << "jfa&egs: propagator is broken\n";
//                (*dit)->update_and_write(bunch);
            }
            ++step_count;
            if (verbose) {
                std::cout << "Propagator:   step " << step_count << "/"
                        << num_steps << std::endl;
            }
            (*it)->apply(bunch);
        }
    }
    for (Multi_diagnostics::iterator it = per_step_diagnostics.begin(); it
            != per_step_diagnostics.end(); ++it) {
        std::cout << "jfa&egs: propagator is broken\n";

        //        (*it)->update_and_write(bunch);
    }
    for (Multi_diagnostics::iterator it = per_turn_diagnostics.begin(); it
            != per_turn_diagnostics.end(); ++it) {
        std::cout << "jfa&egs: propagator is broken\n";
//        (*it)->update_and_write(bunch);
    }
}

Propagator::~Propagator()
{

}
