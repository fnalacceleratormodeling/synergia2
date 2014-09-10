#include "madx_adaptor_map.h"
#include "madx_adaptors.h"

MadX_adaptor_map::MadX_adaptor_map() :
        Element_adaptor_map()
{
    boost::shared_ptr<Marker_madx_adaptor > marker_madx_adaptor(
            new Marker_madx_adaptor);
    set_adaptor("marker", marker_madx_adaptor);
    set_adaptor("endmark", marker_madx_adaptor);

    boost::shared_ptr<Drift_madx_adaptor > drift_madx_adaptor(
            new Drift_madx_adaptor);
    set_adaptor("drift", drift_madx_adaptor);

    boost::shared_ptr<Sbend_madx_adaptor > sbend_madx_adaptor(
            new Sbend_madx_adaptor);
    set_adaptor("sbend", sbend_madx_adaptor);

    boost::shared_ptr<Rbend_madx_adaptor > rbend_madx_adaptor(
            new Rbend_madx_adaptor);
    set_adaptor("rbend", rbend_madx_adaptor);

    boost::shared_ptr<Quadrupole_madx_adaptor > quadrupole_madx_adaptor(
            new Quadrupole_madx_adaptor);
    set_adaptor("quadrupole", quadrupole_madx_adaptor);

    boost::shared_ptr<Sextupole_madx_adaptor > sextupole_madx_adaptor(
            new Sextupole_madx_adaptor);
    set_adaptor("sextupole", sextupole_madx_adaptor);

    boost::shared_ptr<Octupole_madx_adaptor > octupole_madx_adaptor(
            new Octupole_madx_adaptor);
    set_adaptor("octupole", octupole_madx_adaptor);

    boost::shared_ptr<Multipole_madx_adaptor > multipole_madx_adaptor(
            new Multipole_madx_adaptor);
    set_adaptor("multipole", multipole_madx_adaptor);

    boost::shared_ptr<Thinpole_madx_adaptor > thinpole_madx_adaptor(
            new Thinpole_madx_adaptor);
    set_adaptor("thinpole", thinpole_madx_adaptor);

    boost::shared_ptr<Solenoid_madx_adaptor > solenoid_madx_adaptor(
            new Solenoid_madx_adaptor);
    set_adaptor("solenoid", solenoid_madx_adaptor);

    boost::shared_ptr<Hkicker_madx_adaptor > hkicker_madx_adaptor(
            new Hkicker_madx_adaptor);
    set_adaptor("hkicker", hkicker_madx_adaptor);

    boost::shared_ptr<Vkicker_madx_adaptor > vkicker_madx_adaptor(
            new Vkicker_madx_adaptor);
    set_adaptor("vkicker", vkicker_madx_adaptor);

    boost::shared_ptr<Kicker_madx_adaptor > kicker_madx_adaptor(
            new Kicker_madx_adaptor);
    set_adaptor("kicker", kicker_madx_adaptor);

    boost::shared_ptr<Rfcavity_madx_adaptor > rfcavity_madx_adaptor(
            new Rfcavity_madx_adaptor);
    set_adaptor("rfcavity", rfcavity_madx_adaptor);

    boost::shared_ptr<Elseparator_madx_adaptor > elseparator_madx_adaptor(
            new Elseparator_madx_adaptor);
    set_adaptor("elseparator", elseparator_madx_adaptor);

    boost::shared_ptr<Hmonitor_madx_adaptor > hmonitor_madx_adaptor(
            new Hmonitor_madx_adaptor);
    set_adaptor("hmonitor", hmonitor_madx_adaptor);

    boost::shared_ptr<Vmonitor_madx_adaptor > vmonitor_madx_adaptor(
            new Vmonitor_madx_adaptor);
    set_adaptor("vmonitor", vmonitor_madx_adaptor);

    boost::shared_ptr<Monitor_madx_adaptor > monitor_madx_adaptor(
            new Monitor_madx_adaptor);
    set_adaptor("monitor", monitor_madx_adaptor);

    boost::shared_ptr<Instrument_madx_adaptor > instrument_madx_adaptor(
            new Instrument_madx_adaptor);
    set_adaptor("instrument", instrument_madx_adaptor);

    boost::shared_ptr<Ecollimator_madx_adaptor > ecollimator_madx_adaptor(
            new Ecollimator_madx_adaptor);
    set_adaptor("ecollimator", ecollimator_madx_adaptor);

    boost::shared_ptr<Rcollimator_madx_adaptor > rcollimator_madx_adaptor(
            new Rcollimator_madx_adaptor);
    set_adaptor("rcollimator", rcollimator_madx_adaptor);

    boost::shared_ptr<Septum_madx_adaptor > septum_madx_adaptor(
            new Septum_madx_adaptor);
    set_adaptor("e_septum", septum_madx_adaptor);

    boost::shared_ptr<Lambertson_madx_adaptor > lambertson_madx_adaptor(
            new Lambertson_madx_adaptor);
    set_adaptor("lambertson", lambertson_madx_adaptor);
    
     boost::shared_ptr<Srot_madx_adaptor > srot_madx_adaptor(
            new Srot_madx_adaptor);
     set_adaptor("srot", srot_madx_adaptor);

     boost::shared_ptr<Dipedge_madx_adaptor > dipedge_madx_adaptor(
                 new Dipedge_madx_adaptor);
     set_adaptor("dipedge", dipedge_madx_adaptor);
}

template<class Archive>
    void
    MadX_adaptor_map::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element_adaptor_map);
    }

template
void
MadX_adaptor_map::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
MadX_adaptor_map::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
MadX_adaptor_map::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
MadX_adaptor_map::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

BOOST_CLASS_EXPORT_IMPLEMENT(MadX_adaptor_map)
