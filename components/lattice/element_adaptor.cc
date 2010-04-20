#include "element_adaptor.h"

Element_adaptor::Element_adaptor()
{
}

Element_adaptor::~Element_adaptor()
{
}

Element_adaptor_map::Element_adaptor_map()
{
}

void
Element_adaptor_map::set_adaptor(std::string const& name,
        Element_adaptor_sptr const& element_adaptor_sptr)
{
    adaptor_map[name] = element_adaptor_sptr;
}

Element_adaptor_sptr &
Element_adaptor_map::get_adaptor(std::string const& name)
{
    return adaptor_map[name];
}

std::list<std::string >
Element_adaptor_map::get_adaptor_names() const
{
    std::list<std::string > retval;
    for (std::map<std::string, Element_adaptor_sptr >::const_iterator it =
            adaptor_map.begin(); it != adaptor_map.end(); ++it) {
        retval.push_back(it->first);
    }
    return retval;
}

Element_adaptor_map::~Element_adaptor_map()
{
}

Marker_mad8_adaptor::Marker_mad8_adaptor()
{
}

void
Marker_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
}

Marker_mad8_adaptor::~Marker_mad8_adaptor()
{
}

Drift_mad8_adaptor::Drift_mad8_adaptor()
{
}

void
Drift_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
}

Drift_mad8_adaptor::~Drift_mad8_adaptor()
{
}

Sbend_mad8_adaptor::Sbend_mad8_adaptor()
{
}

void
Sbend_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("angle", 0.0);
    lattice_element.set_double_attribute("k1", 0.0);
    lattice_element.set_double_attribute("e1", 0.0);
    lattice_element.set_double_attribute("e2", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
    lattice_element.set_double_attribute("k2", 0.0);
    lattice_element.set_double_attribute("h1", 0.0);
    lattice_element.set_double_attribute("h2", 0.0);
    lattice_element.set_double_attribute("hgap", 0.0);
    lattice_element.set_double_attribute("fint", 0.0);
    lattice_element.set_double_attribute("k3", 0.0);

}

Sbend_mad8_adaptor::~Sbend_mad8_adaptor()
{
}

Rbend_mad8_adaptor::Rbend_mad8_adaptor()
{
}

void
Rbend_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("angle", 0.0);
    lattice_element.set_double_attribute("k1", 0.0);
    lattice_element.set_double_attribute("e1", 0.0);
    lattice_element.set_double_attribute("e2", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
    lattice_element.set_double_attribute("k2", 0.0);
    lattice_element.set_double_attribute("h1", 0.0);
    lattice_element.set_double_attribute("h2", 0.0);
    lattice_element.set_double_attribute("hgap", 0.0);
    lattice_element.set_double_attribute("fint", 0.0);
    lattice_element.set_double_attribute("k3", 0.0);
}

Rbend_mad8_adaptor::~Rbend_mad8_adaptor()
{
}

Quadrupole_mad8_adaptor::Quadrupole_mad8_adaptor()
{
}

void
Quadrupole_mad8_adaptor::set_default_atributes(
        Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("k1", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

Quadrupole_mad8_adaptor::~Quadrupole_mad8_adaptor()
{
}

Sextupole_mad8_adaptor::Sextupole_mad8_adaptor()
{
}

void
Sextupole_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("k2", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

Sextupole_mad8_adaptor::~Sextupole_mad8_adaptor()
{
}

Octupole_mad8_adaptor::Octupole_mad8_adaptor()
{
}

void
Octupole_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("k3", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

Octupole_mad8_adaptor::~Octupole_mad8_adaptor()
{
}

Multipole_mad8_adaptor::Multipole_mad8_adaptor()
{
}

void
Multipole_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("lrad", 0.0);
    lattice_element.set_double_attribute("k0l", 0.0);
    lattice_element.set_double_attribute("t0", 0.0);
    lattice_element.set_double_attribute("k1l", 0.0);
    lattice_element.set_double_attribute("t1", 0.0);
    lattice_element.set_double_attribute("k2l", 0.0);
    lattice_element.set_double_attribute("t2", 0.0);
    lattice_element.set_double_attribute("k3l", 0.0);
    lattice_element.set_double_attribute("t3", 0.0);
    lattice_element.set_double_attribute("k4l", 0.0);
    lattice_element.set_double_attribute("t4", 0.0);
    lattice_element.set_double_attribute("k5l", 0.0);
    lattice_element.set_double_attribute("t5", 0.0);
    lattice_element.set_double_attribute("k6l", 0.0);
    lattice_element.set_double_attribute("t6", 0.0);
    lattice_element.set_double_attribute("k7l", 0.0);
    lattice_element.set_double_attribute("t7", 0.0);
    lattice_element.set_double_attribute("k8l", 0.0);
    lattice_element.set_double_attribute("t8", 0.0);
    lattice_element.set_double_attribute("k9l", 0.0);
    lattice_element.set_double_attribute("t9", 0.0);
}

Multipole_mad8_adaptor::~Multipole_mad8_adaptor()
{
}

Solenoid_mad8_adaptor::Solenoid_mad8_adaptor()
{
}

void
Solenoid_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("ks", 0.0);
}

Solenoid_mad8_adaptor::~Solenoid_mad8_adaptor()
{
}

Hkicker_mad8_adaptor::Hkicker_mad8_adaptor()
{
}

void
Hkicker_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("kick", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

Hkicker_mad8_adaptor::~Hkicker_mad8_adaptor()
{
}

Vkicker_mad8_adaptor::Vkicker_mad8_adaptor()
{
}

void
Vkicker_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("kick", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

Vkicker_mad8_adaptor::~Vkicker_mad8_adaptor()
{
}

Kicker_mad8_adaptor::Kicker_mad8_adaptor()
{
}

void
Kicker_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("hkick", 0.0);
    lattice_element.set_double_attribute("vkick", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

Kicker_mad8_adaptor::~Kicker_mad8_adaptor()
{
}

Rfcavity_mad8_adaptor::Rfcavity_mad8_adaptor()
{
}

void
Rfcavity_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("volt", 0.0);
    lattice_element.set_double_attribute("lag", 0.0);
    lattice_element.set_double_attribute("harmon", 0.0);
    lattice_element.set_double_attribute("betrf", 0.0);
    lattice_element.set_double_attribute("pg", 0.0);
    lattice_element.set_double_attribute("shunt", 0.0);
    lattice_element.set_double_attribute("tfill", 0.0);
}

Rfcavity_mad8_adaptor::~Rfcavity_mad8_adaptor()
{
}

Elseparator_mad8_adaptor::Elseparator_mad8_adaptor()
{
}

void
Elseparator_mad8_adaptor::set_default_atributes(
        Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("e", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

Elseparator_mad8_adaptor::~Elseparator_mad8_adaptor()
{
}

Hmonitor_mad8_adaptor::Hmonitor_mad8_adaptor()
{
}

void
Hmonitor_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

Hmonitor_mad8_adaptor::~Hmonitor_mad8_adaptor()
{
}

Vmonitor_mad8_adaptor::Vmonitor_mad8_adaptor()
{
}

void
Vmonitor_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

Vmonitor_mad8_adaptor::~Vmonitor_mad8_adaptor()
{
}

Monitor_mad8_adaptor::Monitor_mad8_adaptor()
{
}

void
Monitor_mad8_adaptor::set_default_atributes(Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

Monitor_mad8_adaptor::~Monitor_mad8_adaptor()
{
}

Instrument_mad8_adaptor::Instrument_mad8_adaptor()
{
}

void
Instrument_mad8_adaptor::set_default_atributes(
        Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

Instrument_mad8_adaptor::~Instrument_mad8_adaptor()
{
}

Ecollimator_mad8_adaptor::Ecollimator_mad8_adaptor()
{
}

void
Ecollimator_mad8_adaptor::set_default_atributes(
        Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("xsize", 0.0);
    lattice_element.set_double_attribute("ysize", 0.0);
}

Ecollimator_mad8_adaptor::~Ecollimator_mad8_adaptor()
{
}

Rcollimator_mad8_adaptor::Rcollimator_mad8_adaptor()
{
}

void
Rcollimator_mad8_adaptor::set_default_atributes(
        Lattice_element & lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("xsize", 0.0);
    lattice_element.set_double_attribute("ysize", 0.0);
}

Rcollimator_mad8_adaptor::~Rcollimator_mad8_adaptor()
{
}

