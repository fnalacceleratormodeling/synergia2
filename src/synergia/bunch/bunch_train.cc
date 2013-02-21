#include "bunch_train.h"

void
Bunch_train::set_bucket_indices()
{
    std::list<int > found_indices;
    for (int i = 0; i < bunches.size(); ++i) {
        if (bunches[i]->get_bucket_index() == 0) {
            bunches[i]->set_bucket_index(i);
        }
        for (std::list<int >::const_iterator it = found_indices.begin();
                it != found_indices.end(); ++it) {
            if (*it == bunches[i]->get_bucket_index()) {
                throw std::runtime_error(
                        "Bunch_train: bunch bucket indices must be either unique or all zero");
            }
        }
        found_indices.push_back(bunches[i]->get_bucket_index());
    }
}

Bunch_train::Bunch_train(Bunches const& bunches, double spacing) :
                bunches(bunches),
                spacings(std::vector<double >(bunches.size() - 1, spacing))
{
    set_bucket_indices();
}

Bunch_train::Bunch_train(Bunches const& bunches,
        std::vector<double > const& spacings) :
                bunches(bunches),
                spacings(spacings)
{
    if (spacings.size() != bunches.size() - 1) {
        throw std::runtime_error(
                "Bunch_train:: spacings must have length (length(bunches)-1)");
    }
    set_bucket_indices();
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
