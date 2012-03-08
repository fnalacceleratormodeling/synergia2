#include "diagnostics_actions.h"
#include <algorithm>

Diagnostics_actions::Diagnostics_actions() :
    have_bunch_sptr_(false)
{
}

// set_bunches is a local function (template)
template<typename T>
    void
    set_bunch_sptr_all(T & t, Bunch_sptr bunch_sptr)
    {
        for (typename T::iterator it = t.begin(); it != t.end(); ++it) {
            it->diagnostics_sptr->set_bunch_sptr(bunch_sptr);
        }

    }
void
Diagnostics_actions::set_bunch_sptr(Bunch_sptr bunch_sptr)
{
    this->bunch_sptr = bunch_sptr;
    have_bunch_sptr_ = true;

    set_bunch_sptr_all(per_turn_periodic, bunch_sptr);
    set_bunch_sptr_all(per_step_periodic, bunch_sptr);
    set_bunch_sptr_all(per_turn_listed, bunch_sptr);
    set_bunch_sptr_all(per_step_periodic_listed, bunch_sptr);
}

bool
Diagnostics_actions::have_bunch_sptr() const
{
    return have_bunch_sptr_;
}

Bunch_sptr
Diagnostics_actions::get_bunch_sptr()
{
    return bunch_sptr;
}

void
Diagnostics_actions::add_per_turn(Diagnostics_sptr diagnostics_sptr, int period)
{
    if (have_bunch_sptr()) {
        diagnostics_sptr->set_bunch_sptr(get_bunch_sptr());
    }
    per_turn_periodic.push_back(Periodic(period, diagnostics_sptr));
}

void
Diagnostics_actions::add_per_turn(Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& turn_numbers)
{
    if (have_bunch_sptr()) {
        diagnostics_sptr->set_bunch_sptr(get_bunch_sptr());
    }
    per_turn_listed.push_back(Listed(turn_numbers, diagnostics_sptr));
}

void
Diagnostics_actions::add_per_step(Diagnostics_sptr diagnostics_sptr, int period)
{
    if (have_bunch_sptr()) {
        diagnostics_sptr->set_bunch_sptr(get_bunch_sptr());
    }
    per_step_periodic.push_back(Periodic(period, diagnostics_sptr));
}
void
Diagnostics_actions::add_per_step(Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& step_numbers, int turn_period)
{
    if (have_bunch_sptr()) {
        diagnostics_sptr->set_bunch_sptr(get_bunch_sptr());
    }
    per_step_periodic_listed.push_back(
            Periodic_listed(turn_period, step_numbers, diagnostics_sptr));
}

void
Diagnostics_actions::update_and_write_periodics(Periodics & periodics, int num)
{
    for (Periodics::iterator it = periodics.begin(); it != periodics.end(); ++it) {
        if (num % it->period == 0) {
            it->diagnostics_sptr->update_and_write();
        }
    }
}

void
Diagnostics_actions::update_and_write_listeds(Listeds & listeds, int num)
{
    for (Listeds::iterator it = listeds.begin(); it != listeds.end(); ++it) {
        Numbers::iterator pos;
        pos = std::find(it->numbers.begin(), it->numbers.end(), num);
        if (pos != it->numbers.end()) {
            it->diagnostics_sptr->update_and_write();
        }

    }
}

void
Diagnostics_actions::update_and_write_periodic_listeds(
        Periodic_listeds & periodic_listeds, int step, int turn)
{
    for (Periodic_listeds::iterator it = periodic_listeds.begin(); it
            != periodic_listeds.end(); ++it) {
        if (turn % it->turn_period == 0) {
            Numbers::iterator pos;
            pos = std::find(it->step_numbers.begin(), it->step_numbers.end(),
                    step);
            if (pos != it->step_numbers.end()) {
                it->diagnostics_sptr->update_and_write();
            }
        }
    }
}

void
Diagnostics_actions::first_action(Stepper & stepper, Bunch & bunch)
{
    update_and_write_periodics(per_turn_periodic, 0);
    update_and_write_periodics(per_step_periodic, 0);
    update_and_write_listeds(per_turn_listed, 0);
    update_and_write_periodic_listeds(per_step_periodic_listed, 0, 0);
}

void
Diagnostics_actions::turn_end_action(Stepper & stepper, Bunch & bunch,
        int turn_num)
{
    update_and_write_periodics(per_turn_periodic, turn_num);
    update_and_write_listeds(per_turn_listed, turn_num);
}

void
Diagnostics_actions::step_end_action(Stepper & stepper, Step & step,
        Bunch & bunch, int turn_num, int step_num)
{
    update_and_write_periodics(per_step_periodic, step_num);
    update_and_write_periodic_listeds(per_step_periodic_listed, step_num,
            turn_num);
}

Diagnostics_actions::~Diagnostics_actions()
{
}

