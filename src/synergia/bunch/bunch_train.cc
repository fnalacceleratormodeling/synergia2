#include "bunch_train.h"

Bunch_train::Bunch_train(Bunches const& bunches,
        double spacing) :
        bunches(bunches), spacings(
                std::vector<double >(bunches.size() - 1, spacing))
{
}

Bunch_train::Bunch_train(Bunches const& bunches,
        std::vector<double > const& spacings) :
        bunches(bunches), spacings(spacings)
{
    if (spacings.size() != bunches.size() - 1) {
        throw std::runtime_error(
                "Bunch_train:: spacings must have length (length(bunches)-1)");
    }
}

Bunch_train::Bunch_train()
{
}

size_t
Bunch_train::get_size() const
{
    return bunches.size();
}

Bunches &
Bunch_train::get_bunches()
{
    return bunches;
}

std::vector<double > &
Bunch_train::get_spacings()
{
    return spacings;
}

template<class Archive>
    void
    Bunch_train::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(bunches);
        ar & BOOST_SERIALIZATION_NVP(spacings);
    }

template
void
Bunch_train::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Bunch_train::~Bunch_train()
{
}
