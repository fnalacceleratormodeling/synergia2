#include "element_adaptor_map.h"

Element_adaptor_map::Element_adaptor_map() :
        adaptor_map()
{
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

Mad8_adaptor_map::Mad8_adaptor_map() :
        Element_adaptor_map()
{
    boost::shared_ptr<Marker_mad8_adaptor > marker_mad8_adaptor(
            new Marker_mad8_adaptor);
    set_adaptor("marker", marker_mad8_adaptor);

    boost::shared_ptr<Drift_mad8_adaptor > drift_mad8_adaptor(
            new Drift_mad8_adaptor);
    set_adaptor("drift", drift_mad8_adaptor);

    boost::shared_ptr<Sbend_mad8_adaptor > sbend_mad8_adaptor(
            new Sbend_mad8_adaptor);
    set_adaptor("sbend", sbend_mad8_adaptor);

    boost::shared_ptr<Rbend_mad8_adaptor > rbend_mad8_adaptor(
            new Rbend_mad8_adaptor);
    set_adaptor("rbend", rbend_mad8_adaptor);

    boost::shared_ptr<Quadrupole_mad8_adaptor > quadrupole_mad8_adaptor(
            new Quadrupole_mad8_adaptor);
    set_adaptor("quadrupole", quadrupole_mad8_adaptor);

    boost::shared_ptr<Sextupole_mad8_adaptor > sextupole_mad8_adaptor(
            new Sextupole_mad8_adaptor);
    set_adaptor("sextupole", sextupole_mad8_adaptor);

    boost::shared_ptr<Octupole_mad8_adaptor > octupole_mad8_adaptor(
            new Octupole_mad8_adaptor);
    set_adaptor("octupole", octupole_mad8_adaptor);

    boost::shared_ptr<Multipole_mad8_adaptor > multipole_mad8_adaptor(
            new Multipole_mad8_adaptor);
    set_adaptor("multipole", multipole_mad8_adaptor);

    boost::shared_ptr<Thinpole_mad8_adaptor > thinpole_mad8_adaptor(
            new Thinpole_mad8_adaptor);
    set_adaptor("thinpole", thinpole_mad8_adaptor);

    boost::shared_ptr<Solenoid_mad8_adaptor > solenoid_mad8_adaptor(
            new Solenoid_mad8_adaptor);
    set_adaptor("solenoid", solenoid_mad8_adaptor);

    boost::shared_ptr<Hkicker_mad8_adaptor > hkicker_mad8_adaptor(
            new Hkicker_mad8_adaptor);
    set_adaptor("hkicker", hkicker_mad8_adaptor);

    boost::shared_ptr<Vkicker_mad8_adaptor > vkicker_mad8_adaptor(
            new Vkicker_mad8_adaptor);
    set_adaptor("vkicker", vkicker_mad8_adaptor);

    boost::shared_ptr<Kicker_mad8_adaptor > kicker_mad8_adaptor(
            new Kicker_mad8_adaptor);
    set_adaptor("kicker", kicker_mad8_adaptor);

    boost::shared_ptr<Rfcavity_mad8_adaptor > rfcavity_mad8_adaptor(
            new Rfcavity_mad8_adaptor);
    set_adaptor("rfcavity", rfcavity_mad8_adaptor);

    boost::shared_ptr<Elseparator_mad8_adaptor > elseparator_mad8_adaptor(
            new Elseparator_mad8_adaptor);
    set_adaptor("elseparator", elseparator_mad8_adaptor);

    boost::shared_ptr<Hmonitor_mad8_adaptor > hmonitor_mad8_adaptor(
            new Hmonitor_mad8_adaptor);
    set_adaptor("hmonitor", hmonitor_mad8_adaptor);

    boost::shared_ptr<Vmonitor_mad8_adaptor > vmonitor_mad8_adaptor(
            new Vmonitor_mad8_adaptor);
    set_adaptor("vmonitor", vmonitor_mad8_adaptor);

    boost::shared_ptr<Monitor_mad8_adaptor > monitor_mad8_adaptor(
            new Monitor_mad8_adaptor);
    set_adaptor("monitor", monitor_mad8_adaptor);

    boost::shared_ptr<Instrument_mad8_adaptor > instrument_mad8_adaptor(
            new Instrument_mad8_adaptor);
    set_adaptor("instrument", instrument_mad8_adaptor);

    boost::shared_ptr<Ecollimator_mad8_adaptor > ecollimator_mad8_adaptor(
            new Ecollimator_mad8_adaptor);
    set_adaptor("ecollimator", ecollimator_mad8_adaptor);

    boost::shared_ptr<Rcollimator_mad8_adaptor > rcollimator_mad8_adaptor(
            new Rcollimator_mad8_adaptor);
    set_adaptor("rcollimator", rcollimator_mad8_adaptor);

    boost::shared_ptr<Septum_mad8_adaptor > septum_mad8_adaptor(
            new Septum_mad8_adaptor);
    set_adaptor("e_septum", septum_mad8_adaptor);

    boost::shared_ptr<Lambertson_mad8_adaptor > lambertson_mad8_adaptor(
            new Lambertson_mad8_adaptor);
    set_adaptor("lambertson", lambertson_mad8_adaptor);
}

template<class Archive>
    void
    Mad8_adaptor_map::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor_map);
    }

template
void
Mad8_adaptor_map::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Mad8_adaptor_map::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Mad8_adaptor_map::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Mad8_adaptor_map::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

BOOST_CLASS_EXPORT_IMPLEMENT(Mad8_adaptor_map)
