#include "ff_element_map.h"

FF_element_map::FF_element_map()
{

}

bool FF_element_map::has_element_type(std::string const& type) const
{
    return (element_map.count(type) > 0);
}

void FF_element_map::set_element_type(std::string const& type,
                             FF_element_sptr element_sptr)
{
    element_map[type] = element_sptr;
}

FF_element_sptr FF_element_map::get_element_type(std::string const& type) const
{
    std::map<std::string, FF_element_sptr >::const_iterator iter =
            element_map.find(type);
    if (iter == element_map.end()) {
        throw std::runtime_error("FF_element_map::get_element_type: type '" + type + "' not found");
    }
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

// jfa "the_big_giant_global_ff_element_map" is a temporary workaround for the
// problem of making the FF_element_map available to FF_propagate_operation.
#include "ff_drift.h"
#include "ff_sbend.h"
#include "ff_rbend.h"
#include "ff_multipole.h"
#include "ff_quadrupole.h"
#include "ff_sextupole.h"
#include "ff_octupole.h"
#include "ff_rfcavity.h"
#include "ff_kicker.h"
#include "ff_hkicker.h"
#include "ff_vkicker.h"
#include "ff_marker.h"
#include "ff_constfoc.h"
#include "ff_nllens.h"
#include "ff_dipedge.h"
#include "ff_solenoid.h"
#include "ff_elens.h"
#include "ff_misc_elements.h"

FF_element_map the_big_giant_global_ff_element_map;
void construct_big_giant_global_ff_element_map()
{
    the_big_giant_global_ff_element_map.set_element_type("drift",      FF_drift_sptr(new FF_drift()));
    the_big_giant_global_ff_element_map.set_element_type("sbend",      FF_sbend_sptr(new FF_sbend()));
    the_big_giant_global_ff_element_map.set_element_type("rbend",      FF_rbend_sptr(new FF_rbend()));
    the_big_giant_global_ff_element_map.set_element_type("rfcavity",   FF_rfcavity_sptr(new FF_rfcavity()));
    the_big_giant_global_ff_element_map.set_element_type("kicker",     FF_kicker_sptr(new FF_kicker()));
    the_big_giant_global_ff_element_map.set_element_type("hkicker",    FF_hkicker_sptr(new FF_hkicker()));
    the_big_giant_global_ff_element_map.set_element_type("vkicker",    FF_vkicker_sptr(new FF_vkicker()));
    the_big_giant_global_ff_element_map.set_element_type("multipole",  FF_multipole_sptr(new FF_multipole()));
    the_big_giant_global_ff_element_map.set_element_type("quadrupole", FF_quadrupole_sptr(new FF_quadrupole()));
    the_big_giant_global_ff_element_map.set_element_type("sextupole",  FF_sextupole_sptr(new FF_sextupole()));
    the_big_giant_global_ff_element_map.set_element_type("octupole",   FF_octupole_sptr(new FF_octupole()));
    the_big_giant_global_ff_element_map.set_element_type("marker",     FF_marker_sptr(new FF_marker()));
    the_big_giant_global_ff_element_map.set_element_type("monitor",    FF_monitor_sptr(new FF_monitor()));
    the_big_giant_global_ff_element_map.set_element_type("hmonitor",   FF_hmonitor_sptr(new FF_hmonitor()));
    the_big_giant_global_ff_element_map.set_element_type("vmonitor",   FF_vmonitor_sptr(new FF_vmonitor()));
    the_big_giant_global_ff_element_map.set_element_type("constfoc",   FF_constfoc_sptr(new FF_constfoc()));
    the_big_giant_global_ff_element_map.set_element_type("nllens",     FF_nllens_sptr(new FF_nllens()));
    the_big_giant_global_ff_element_map.set_element_type("dipedge",    FF_dipedge_sptr(new FF_dipedge()));
    the_big_giant_global_ff_element_map.set_element_type("solenoid",   FF_solenoid_sptr(new FF_solenoid()));
    the_big_giant_global_ff_element_map.set_element_type("elens",      FF_elens_sptr(new FF_elens()));
    the_big_giant_global_ff_element_map.set_element_type("instrument", FF_instrument_sptr(new FF_instrument()));
}
