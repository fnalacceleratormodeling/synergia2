#include "bunch_simulator.h"
#include "diagnostics_actions.h"

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr) :
    bunch_sptr(bunch_sptr),
            diagnostics_actions_sptr(new Diagnostics_actions)
{
    diagnostics_actions_sptr->set_bunch_sptr(bunch_sptr);
}

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr,
        Diagnostics_actions_sptr diagnostics_actions_sptr) :
    bunch_sptr(bunch_sptr), diagnostics_actions_sptr(diagnostics_actions_sptr)
{
    diagnostics_actions_sptr->set_bunch_sptr(bunch_sptr);
}

Bunch_simulator::Bunch_simulator()
{
}

Bunch &
Bunch_simulator::get_bunch()
{
    return *bunch_sptr;
}

Bunch_sptr
Bunch_simulator::get_bunch_sptr()
{
    return bunch_sptr;
}

Diagnostics_actions &
Bunch_simulator::get_diagnostics_actions()
{
    return *diagnostics_actions_sptr;
}

Diagnostics_actions_sptr
Bunch_simulator::get_diagnostics_actions_sptr()
{
    return diagnostics_actions_sptr;
}

void
Bunch_simulator::add_per_turn(Diagnostics_sptr diagnostics_sptr, int period)
{
    diagnostics_actions_sptr->add_per_turn(diagnostics_sptr, period);
}

void
Bunch_simulator::add_per_turn(Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& turn_numbers)
{
    diagnostics_actions_sptr->add_per_turn(diagnostics_sptr, turn_numbers);
}

void
Bunch_simulator::add_per_step(Diagnostics_sptr diagnostics_sptr, int period)
{
    diagnostics_actions_sptr->add_per_step(diagnostics_sptr, period);
}

void
Bunch_simulator::add_per_step(Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& step_numbers, int turn_period)
{
    diagnostics_actions_sptr->add_per_step(diagnostics_sptr, step_numbers,
            turn_period);
}

void
Bunch_simulator::add_per_forced_diagnostics_step(
        Diagnostics_sptr diagnostics_sptr, int turn_period)
{
    diagnostics_actions_sptr->add_per_forced_diagnostics_step(diagnostics_sptr,
            turn_period);
}

template<class Archive>
    void
    Bunch_simulator::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(bunch_sptr);
        ar & BOOST_SERIALIZATION_NVP(diagnostics_actions_sptr);
    }

template
void
Bunch_simulator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Bunch_simulator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Bunch_simulator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Bunch_simulator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Bunch_simulator::~Bunch_simulator()
{
}
