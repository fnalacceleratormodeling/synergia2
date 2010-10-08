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

void
Propagator::propagate(Bunch & bunch, int num_turns,
        Diagnostics_writer & per_step_diagnostics,
        Diagnostics_writer & per_turn_diagnostics)
{
    for (int turn = 0; turn < num_turns; ++turn) {
        std::cout << "Propagator: turn " << turn+1 << "/" << num_turns
                << std::endl;
        bunch.get_reference_particle().start_repetition();
        per_turn_diagnostics.update_and_write(bunch);
        int step_count = 0;
        for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
                != stepper_sptr->get_steps().end(); ++it) {
            per_step_diagnostics.update_and_write(bunch);
            ++step_count;
            std::cout << "Propagator: step " << step_count << std::endl;
            (*it)->apply(bunch);
        }
    }
    per_step_diagnostics.update_and_write(bunch);
    per_turn_diagnostics.update_and_write(bunch);
}

Propagator::~Propagator()
{

}
