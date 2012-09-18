#include "bunch_train_simulator.h"

Bunch_train_simulator::Bunch_train_simulator(Bunches const& bunches,
        double spacing) :
        bunches(bunches), spacings(
                std::vector<double >(bunches.size() - 1, spacing)), diagnostics_actionss(
                bunches.size())
{
}

Bunch_train_simulator::Bunch_train_simulator(Bunches const& bunches,
        std::vector<double > const& spacings) :
        bunches(bunches), spacings(spacings), diagnostics_actionss(
                bunches.size())
{
    if (spacings.size() != bunches.size() - 1) {
        throw std::runtime_error(
                "Bunch_train_simulator:: spacings must have length (length(bunches)-1)");
    }
}

Bunch_train_simulator::Bunch_train_simulator()
{
}

size_t
Bunch_train_simulator::get_size() const
{
    return bunches.size();
}

Bunches &
Bunch_train_simulator::get_bunches()
{
    return bunches;
}


std::vector<double > &
Bunch_train_simulator::get_spacings()
{
    return spacings;
}

Diagnostics_actionss &
Bunch_train_simulator::get_diagnostics_actionss()
{
    return diagnostics_actionss;
}

void
Bunch_train_simulator::add_per_turn(int which, Diagnostics_sptr diagnostics_sptr,
        int period)
{
    diagnostics_actionss.at(which)->add_per_turn(diagnostics_sptr, period);
}

void
Bunch_train_simulator::add_per_turn(int which, Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& turn_numbers)
{
    diagnostics_actionss.at(which)->add_per_turn(diagnostics_sptr, turn_numbers);
}

void
Bunch_train_simulator::add_per_step(int which, Diagnostics_sptr diagnostics_sptr,
        int period)
{
    diagnostics_actionss.at(which)->add_per_step(diagnostics_sptr, period);
}

void
Bunch_train_simulator::add_per_step(int which, Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& step_numbers, int turn_period)
{
    diagnostics_actionss.at(which)->add_per_step(diagnostics_sptr, step_numbers,
            turn_period);
}

void
Bunch_train_simulator::add_per_forced_diagnostics_step(int which,
        Diagnostics_sptr diagnostics_sptr, int turn_period)
{
    diagnostics_actionss.at(which)->add_per_forced_diagnostics_step(diagnostics_sptr,
            turn_period);
}

template<class Archive>
    void
    Bunch_train_simulator::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(bunches);
        ar & BOOST_SERIALIZATION_NVP(spacings);
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
