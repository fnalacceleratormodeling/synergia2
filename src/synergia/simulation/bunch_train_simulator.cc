#include "bunch_train_simulator.h"

Bunch_train_simulator::Bunch_train_simulator(Bunch_train_sptr bunch_train_sptr) :
        bunch_train_sptr(bunch_train_sptr), diagnostics_actionss(
                bunch_train_sptr->get_size())
{
    for (int i = 0; i < bunch_train_sptr->get_size(); ++i) {
        diagnostics_actionss.push_back(
                Diagnostics_actions_sptr(new Diagnostics_actions));
        diagnostics_actionss.at(i)->set_bunch_sptr(
                bunch_train_sptr->get_bunches().at(i));
    }
}

Bunch_train_simulator::Bunch_train_simulator()
{
}

Bunch_train &
Bunch_train_simulator::get_bunch_train()
{
    return *bunch_train_sptr;
}

Bunch_train_sptr
Bunch_train_simulator::get_bunch_train_sptr()
{
    return bunch_train_sptr;
}

Diagnostics_actionss &
Bunch_train_simulator::get_diagnostics_actionss()
{
    return diagnostics_actionss;
}

void
Bunch_train_simulator::add_per_turn(int which,
        Diagnostics_sptr diagnostics_sptr, int period)
{
    diagnostics_actionss.at(which)->add_per_turn(diagnostics_sptr, period);
}

void
Bunch_train_simulator::add_per_turn(int which,
        Diagnostics_sptr diagnostics_sptr, std::list<int > const& turn_numbers)
{
    diagnostics_actionss.at(which)->add_per_turn(diagnostics_sptr,
            turn_numbers);
}

void
Bunch_train_simulator::add_per_step(int which,
        Diagnostics_sptr diagnostics_sptr, int period)
{
    diagnostics_actionss.at(which)->add_per_step(diagnostics_sptr, period);
}

void
Bunch_train_simulator::add_per_step(int which,
        Diagnostics_sptr diagnostics_sptr, std::list<int > const& step_numbers,
        int turn_period)
{
    diagnostics_actionss.at(which)->add_per_step(diagnostics_sptr, step_numbers,
            turn_period);
}

void
Bunch_train_simulator::add_per_forced_diagnostics_step(int which,
        Diagnostics_sptr diagnostics_sptr, int turn_period)
{
    diagnostics_actionss.at(which)->add_per_forced_diagnostics_step(
            diagnostics_sptr, turn_period);
}

template<class Archive>
    void
    Bunch_train_simulator::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(bunch_train_sptr);
        ar & BOOST_SERIALIZATION_NVP(diagnostics_actionss);
    }

template
void
Bunch_train_simulator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Bunch_train_simulator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Bunch_train_simulator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Bunch_train_simulator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Bunch_train_simulator::~Bunch_train_simulator()
{
}
