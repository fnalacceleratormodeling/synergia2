#include "diagnostics_actions.h"
#include "synergia/utils/string_utils.h"
#include <algorithm>

template<class Archive>
    void
    Diagnostics_actions::Periodic::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(period);
        ar & BOOST_SERIALIZATION_NVP(diagnostics_sptr);
    }

template
void
Diagnostics_actions::Periodic::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Periodic::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Periodic::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Periodic::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_actions::Diagnostics_actions() :
    have_bunch_sptr_(false)
{
}

template<class Archive>
    void
    Diagnostics_actions::Listed::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(numbers);
        ar & BOOST_SERIALIZATION_NVP(diagnostics_sptr);
    }

template
void
Diagnostics_actions::Listed::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Listed::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Listed::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Listed::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

template<class Archive>
    void
    Diagnostics_actions::Periodic_listed::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(turn_period);
        ar & BOOST_SERIALIZATION_NVP(step_numbers);
        ar & BOOST_SERIALIZATION_NVP(diagnostics_sptr);
    }

template
void
Diagnostics_actions::Periodic_listed::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Periodic_listed::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Periodic_listed::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::Periodic_listed::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

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
    set_bunch_sptr_all(per_forced_step_periodic, bunch_sptr);
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
Diagnostics_actions::add_per_turn(Diagnostics_sptr diagnostics_sptr,
        int turn_period)
{
    if (have_bunch_sptr()) {
        diagnostics_sptr->set_bunch_sptr(get_bunch_sptr());
    }
    per_turn_periodic.push_back(Periodic(turn_period, diagnostics_sptr));
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
Diagnostics_actions::add_per_step(Diagnostics_sptr diagnostics_sptr,
        int step_period)
{
    if (have_bunch_sptr()) {
        diagnostics_sptr->set_bunch_sptr(get_bunch_sptr());
    }
    per_step_periodic.push_back(Periodic(step_period, diagnostics_sptr));
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
Diagnostics_actions::add_per_forced_diagnostics_step(
        Diagnostics_sptr diagnostics_sptr, int turn_period)
{
    if (have_bunch_sptr()) {
        diagnostics_sptr->set_bunch_sptr(get_bunch_sptr());
    }
    per_forced_step_periodic.push_back(Periodic(turn_period, diagnostics_sptr));
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
    Listeds::iterator it = listeds.begin();
    while (it != listeds.end()) {
        // what is the maximum turn in this Listed instance?
        int maxturn = -1;
        for (Numbers::iterator tlit=it->numbers.begin();
             tlit!=it->numbers.end(); ++tlit) {
            if (*tlit > maxturn) maxturn = *tlit;
        }
        // remove this particular Listed instance if we have already gone
        // past its maximum turn
        if (num > maxturn) {
            it->diagnostics_sptr->delete_write_helper_ptr(); // close file
            listeds.erase(it++); // erase on a list  preserves existing iterators
        } else {
            Numbers::iterator pos;
            pos = std::find(it->numbers.begin(), it->numbers.end(), num);
            if (pos != it->numbers.end()) {
                it->diagnostics_sptr->update_and_write();
            }
            ++it;
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

    if (per_forced_step_periodic.size() > 0) {
        bool force_diagnostics = false;
        for (Operators::iterator oit = step.get_operators().begin();
                oit != step.get_operators().end(); ++oit) {
            if ((*oit)->get_type() == "independent") {
                Lattice_element_slices slices(
                        boost::static_pointer_cast<Independent_operator >(*oit)->get_slices());
                for (Lattice_element_slices::iterator slit = slices.begin();
                        slit != slices.end(); ++slit) {
                    if ((*slit)->has_right_edge()
                            && (*slit)->get_lattice_element().has_string_attribute(
                                    Stepper::force_diagnostics_attribute)) {
                        if (!false_string(
                                (*slit)->get_lattice_element().get_string_attribute(
                                        Stepper::force_diagnostics_attribute))) {
                            force_diagnostics = true;
                        }
                    }
                }
            }
        }
        if (force_diagnostics) {
            update_and_write_periodics(per_forced_step_periodic, turn_num);
        }
    }
}

template<class Archive>
    void
    Diagnostics_actions::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(bunch_sptr);
        ar & BOOST_SERIALIZATION_NVP(have_bunch_sptr_);
        ar & BOOST_SERIALIZATION_NVP(per_turn_periodic);
        ar & BOOST_SERIALIZATION_NVP(per_step_periodic);
        ar & BOOST_SERIALIZATION_NVP(per_forced_step_periodic);
        ar & BOOST_SERIALIZATION_NVP(per_turn_listed);
        ar & BOOST_SERIALIZATION_NVP(per_step_periodic_listed);
    }

template
void
Diagnostics_actions::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_actions::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_actions::~Diagnostics_actions()
{
}

