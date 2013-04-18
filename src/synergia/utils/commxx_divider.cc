#include "commxx_divider.h"
#include "containers_to_string.h"
#include <stdexcept>

Commxx_divider::Commxx_divider() :
        cache(), subsize(0), per_host(false)
{
}

Commxx_divider::Commxx_divider(int subsize, bool per_host) :
        cache(), subsize(subsize), per_host(per_host)
{
}

Commxx_sptr
Commxx_divider::get_commxx_sptr(Commxx_sptr const& parent)
{
    Commxx_sptr retval;
    std::map<Commxx_sptr, Commxx_sptr >::iterator pos(cache.find(parent));
    if (pos == cache.end()) {
        int parent_size = parent->get_size();
        if ((subsize == 0) || (subsize >= parent_size)) {
            retval = Commxx_sptr(new Commxx(parent, per_host));
            cache[parent] = retval;
        } else {
            if (parent_size % subsize != 0) {
                throw std::runtime_error(
                        "Commxx_divider: parent communicator size must be an integer multiple of subsize");
            }
            std::vector<int > ranks(subsize);
            int min_rank = (parent->get_rank() / subsize) * subsize;
            for (int i = 0; i < subsize; ++i) {
                ranks[i] = min_rank + i;
            }
            retval = Commxx_sptr(new Commxx(parent, ranks, per_host));
            cache[parent] = retval;
        }
    } else {
        retval = pos->second;
    }
    return retval;
}

template<class Archive>
    void
    Commxx_divider::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(cache);
    ar & BOOST_SERIALIZATION_NVP(subsize);
    ar & BOOST_SERIALIZATION_NVP(per_host);
}

template
void
Commxx_divider::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Commxx_divider::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Commxx_divider::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Commxx_divider::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);


Commxx_divider::~Commxx_divider()
{
}
