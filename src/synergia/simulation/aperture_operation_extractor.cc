#include "aperture_operation_extractor.h"

BOOST_CLASS_EXPORT_IMPLEMENT(Circular_extractor)
BOOST_CLASS_EXPORT_IMPLEMENT(Elliptical_extractor)
BOOST_CLASS_EXPORT_IMPLEMENT(Rectangular_extractor)
BOOST_CLASS_EXPORT_IMPLEMENT(Polygon_extractor)
BOOST_CLASS_EXPORT_IMPLEMENT(Wire_elliptical_extractor)
BOOST_CLASS_EXPORT_IMPLEMENT(Lambertson_extractor)
BOOST_CLASS_EXPORT_IMPLEMENT(Rectangular_with_ears_extractor)

Aperture_operation_extractor::Aperture_operation_extractor()
{
}

template <class Archive>
void Aperture_operation_extractor::serialize(Archive & ar, const unsigned int version) {
}

template
void
Aperture_operation_extractor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Aperture_operation_extractor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Aperture_operation_extractor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Aperture_operation_extractor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Aperture_operation_extractor::~Aperture_operation_extractor()
{
}

Aperture_operation_extractor_map::Aperture_operation_extractor_map()
{
}

void
Aperture_operation_extractor_map::set_extractor(std::string const& name,
        Aperture_operation_extractor_sptr operation_extractor)
{
    extractor_map[name] = operation_extractor;
}

Aperture_operation_extractor_sptr
Aperture_operation_extractor_map::get_extractor(std::string const& name)
{
    return extractor_map[name];
}

std::list<std::string >
Aperture_operation_extractor_map::get_extractor_names() const
{
    std::list<std::string > retval;
    for (std::map<std::string, Aperture_operation_extractor_sptr >::const_iterator
            it = extractor_map.begin(); it != extractor_map.end(); ++it) {
        retval.push_back(it->first);
    }
    return retval;
}

template <class Archive>
void Aperture_operation_extractor_map::serialize(Archive & ar, const unsigned int version) {
    ar & BOOST_SERIALIZATION_NVP(extractor_map);
}

template
void
Aperture_operation_extractor_map::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Aperture_operation_extractor_map::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Aperture_operation_extractor_map::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Aperture_operation_extractor_map::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Aperture_operation_extractor_map::~Aperture_operation_extractor_map()
{
}
