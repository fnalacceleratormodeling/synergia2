#include "bunch_with_diagnostics.h"
#include "synergia/simulation/standard_diagnostics_actions.cc"
#include "synergia/simulation/propagate_actions.cc"
#include <boost/pointer_cast.hpp>

Bunch_with_diagnostics::Bunch_with_diagnostics(Bunch_sptr bunch_sptr,
        Standard_diagnostics_actions_sptr diagnostics_actions_sptr) :
    bunch_sptr(bunch_sptr), diagnostics_actions_sptr(diagnostics_actions_sptr)
{
}

Bunch_with_diagnostics::Bunch_with_diagnostics()
{
}

void
Bunch_with_diagnostics::check_bunch_pointer_in_diagnostics() const
{

    const std::list<Diagnostics_sptr > steps_list =
            get_diagnostics_actions_sptr()->get_per_steps_diagnostics_list();
    if (!steps_list.empty()) {
        for (std::list<Diagnostics_sptr >::const_iterator it =
                steps_list.begin(); it != steps_list.end(); ++it) {
            Diagnostics_sptr diagnostics_sptr(
                    boost::dynamic_pointer_cast<Diagnostics >(*it));
            if (diagnostics_sptr->get_bunch_sptr() != bunch_sptr) throw std::runtime_error(
                    "bunch_with_diagnostics: the step_diagnostics bunch pointer and the bunch are different objects");

        }

    }

    const std::list<Diagnostics_sptr > turns_list =
            get_diagnostics_actions_sptr()->get_per_turns_diagnostics_list();
    if (!turns_list.empty()) {
        for (std::list<Diagnostics_sptr >::const_iterator it =
                turns_list.begin(); it != turns_list.end(); ++it) {
            Diagnostics_sptr diagnostics_sptr(
                    boost::dynamic_pointer_cast<Diagnostics >(*it));
            if (diagnostics_sptr->get_bunch_sptr() != bunch_sptr) throw std::runtime_error(
                    "bunch_with_diagnostics: the turn_diagnostics bunch pointer and the bunch are different objects");
        }
    }

}

void
Bunch_with_diagnostics::add_per_step_diagnostics(
        Diagnostics_sptr diagnostics_sptr)
{
    get_diagnostics_actions_sptr()->add_per_step(diagnostics_sptr);
    check_bunch_pointer_in_diagnostics();
}

void
Bunch_with_diagnostics::add_per_turn_diagnostics(
        Diagnostics_sptr diagnostics_sptr)
{
    get_diagnostics_actions_sptr()->add_per_turn(diagnostics_sptr);
    check_bunch_pointer_in_diagnostics();
}

Bunch_sptr const
Bunch_with_diagnostics::get_bunch_sptr()
{
    return bunch_sptr;
}

Standard_diagnostics_actions_sptr
Bunch_with_diagnostics::get_diagnostics_actions_sptr() const
{
    return diagnostics_actions_sptr;
}

Commxx const&
Bunch_with_diagnostics::get_comm() const
{
    return bunch_sptr->get_comm();
}

Bunch_with_diagnostics::~Bunch_with_diagnostics()
{
}

/* */
