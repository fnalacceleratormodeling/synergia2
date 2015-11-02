#include "element_adaptor_map.h"

Element_adaptor_map::Element_adaptor_map() :
        adaptor_map()
      , label("")
{
}

Element_adaptor_map::Element_adaptor_map(std::string const& label) :
    adaptor_map()
  , label(label)
{
}

std::string
Element_adaptor_map::get_label() const
{
    return label;
}

void
Element_adaptor_map::set_adaptor(std::string const& name,
        Element_adaptor_sptr element_adaptor_sptr)
{
    adaptor_map[name] = element_adaptor_sptr;
}

bool
Element_adaptor_map::has_adaptor(std::string const& name) const
{
    return (adaptor_map.count(name) > 0);
}

Element_adaptor_sptr
Element_adaptor_map::get_adaptor(std::string const& name) const
{
    std::map<std::string, Element_adaptor_sptr >::const_iterator iter =
            adaptor_map.find(name);
    return iter->second;
}

std::list<std::string >
Element_adaptor_map::get_adaptor_names() const
{
    std::list < std::string > retval;
    for (std::map<std::string, Element_adaptor_sptr >::const_iterator it =
            adaptor_map.begin(); it != adaptor_map.end(); ++it) {
        retval.push_back(it->first);
    }
    return retval;
}

template<class Archive>
    void
    Element_adaptor_map::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(adaptor_map);
        ar & BOOST_SERIALIZATION_NVP(label);
    }

template
void
Element_adaptor_map::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Element_adaptor_map::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Element_adaptor_map::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Element_adaptor_map::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Element_adaptor_map::~Element_adaptor_map()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Element_adaptor_map)
