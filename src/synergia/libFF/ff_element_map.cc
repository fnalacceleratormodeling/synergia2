#include "ff_element_map.h"

FF_element_map::FF_element_map()
{

}

bool FF_element_map::has_adaptor(std::string const& type) const
{
    return (element_map.count(type) > 0);
}

void FF_element_map::set_adaptor(std::string const& type,
                             FF_element_sptr element_sptr)
{
    element_map[type] = element_sptr;
}

FF_element_sptr FF_element_map::get_adaptor(std::string const& type) const
{
    std::map<std::string, FF_element_sptr >::const_iterator iter =
            element_map.find(type);
    return iter->second;
}

template<class Archive>
    void
    FF_element_map::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(element_map);
    }

template
void
FF_element_map::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_element_map::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_element_map::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_element_map::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_element_map::~FF_element_map()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(FF_element_map)
