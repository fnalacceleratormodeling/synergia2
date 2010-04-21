#include "element_adaptor.h"

Element_adaptor::Element_adaptor()
{
}

void
Element_adaptor::set_double_default(Lattice_element & lattice_element,
        std::string const& name, double value)
{
    if (!lattice_element.has_double_attribute(name)) {
        lattice_element.set_double_attribute(name, value);
    }
}

void
set_string_default(Lattice_element & lattice_element, std::string const& name,
        std::string const& value)
{
    if (!lattice_element.has_string_attribute(name)) {
        lattice_element.set_string_attribute(name, value);
    }
}

Element_adaptor::~Element_adaptor()
{
}

Element_adaptor_map::Element_adaptor_map()
{
    boost::shared_ptr<Marker_mad8_adaptor > marker_mad8_adaptor(
            new Marker_mad8_adaptor);
    adaptor_map["marker"] = marker_mad8_adaptor;

    boost::shared_ptr<Drift_mad8_adaptor > drift_mad8_adaptor(
            new Drift_mad8_adaptor);
    adaptor_map["drift"] = drift_mad8_adaptor;

    boost::shared_ptr<Sbend_mad8_adaptor > sbend_mad8_adaptor(
            new Sbend_mad8_adaptor);
    adaptor_map["sbend"] = sbend_mad8_adaptor;

    boost::shared_ptr<Rbend_mad8_adaptor > rbend_mad8_adaptor(
            new Rbend_mad8_adaptor);
    adaptor_map["rbend"] = rbend_mad8_adaptor;

    boost::shared_ptr<Quadrupole_mad8_adaptor > quadrupole_mad8_adaptor(
            new Quadrupole_mad8_adaptor);
    adaptor_map["quadrupole"] = quadrupole_mad8_adaptor;

    boost::shared_ptr<Sextupole_mad8_adaptor > sextupole_mad8_adaptor(
            new Sextupole_mad8_adaptor);
    adaptor_map["sextupole"] = sextupole_mad8_adaptor;

    boost::shared_ptr<Octupole_mad8_adaptor > octupole_mad8_adaptor(
            new Octupole_mad8_adaptor);
    adaptor_map["octupole"] = octupole_mad8_adaptor;

    boost::shared_ptr<Multipole_mad8_adaptor > multipole_mad8_adaptor(
            new Multipole_mad8_adaptor);
    adaptor_map["multipole"] = multipole_mad8_adaptor;

    boost::shared_ptr<Solenoid_mad8_adaptor > solenoid_mad8_adaptor(
            new Solenoid_mad8_adaptor);
    adaptor_map["solenoid"] = solenoid_mad8_adaptor;

    boost::shared_ptr<Hkicker_mad8_adaptor > hkicker_mad8_adaptor(
            new Hkicker_mad8_adaptor);
    adaptor_map["hkicker"] = hkicker_mad8_adaptor;

    boost::shared_ptr<Vkicker_mad8_adaptor > vkicker_mad8_adaptor(
            new Vkicker_mad8_adaptor);
    adaptor_map["vkicker"] = vkicker_mad8_adaptor;

    boost::shared_ptr<Kicker_mad8_adaptor > kicker_mad8_adaptor(
            new Kicker_mad8_adaptor);
    adaptor_map["kicker"] = kicker_mad8_adaptor;

    boost::shared_ptr<Rfcavity_mad8_adaptor > rfcavity_mad8_adaptor(
            new Rfcavity_mad8_adaptor);
    adaptor_map["rfcavity"] = rfcavity_mad8_adaptor;

    boost::shared_ptr<Elseparator_mad8_adaptor > elseparator_mad8_adaptor(
            new Elseparator_mad8_adaptor);
    adaptor_map["elseparator"] = elseparator_mad8_adaptor;

    boost::shared_ptr<Hmonitor_mad8_adaptor > hmonitor_mad8_adaptor(
            new Hmonitor_mad8_adaptor);
    adaptor_map["hmonitor"] = hmonitor_mad8_adaptor;

    boost::shared_ptr<Vmonitor_mad8_adaptor > vmonitor_mad8_adaptor(
            new Vmonitor_mad8_adaptor);
    adaptor_map["vmonitor"] = vmonitor_mad8_adaptor;

    boost::shared_ptr<Monitor_mad8_adaptor > monitor_mad8_adaptor(
            new Monitor_mad8_adaptor);
    adaptor_map["monitor"] = monitor_mad8_adaptor;

    boost::shared_ptr<Instrument_mad8_adaptor > instrument_mad8_adaptor(
            new Instrument_mad8_adaptor);
    adaptor_map["instrument"] = instrument_mad8_adaptor;

    boost::shared_ptr<Ecollimator_mad8_adaptor > ecollimator_mad8_adaptor(
            new Ecollimator_mad8_adaptor);
    adaptor_map["ecollimator"] = ecollimator_mad8_adaptor;

    boost::shared_ptr<Rcollimator_mad8_adaptor > rcollimator_mad8_adaptor(
            new Rcollimator_mad8_adaptor);
    adaptor_map["rcollimator"] = rcollimator_mad8_adaptor;
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
Marker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
}

Marker_mad8_adaptor::~Marker_mad8_adaptor()
{
}

Drift_mad8_adaptor::Drift_mad8_adaptor()
{
}

void
Drift_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
}

Drift_mad8_adaptor::~Drift_mad8_adaptor()
{
}

Sbend_mad8_adaptor::Sbend_mad8_adaptor()
{
}

void
Sbend_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "angle", 0.0);
    set_double_default(lattice_element, "k1", 0.0);
    set_double_default(lattice_element, "e1", 0.0);
    set_double_default(lattice_element, "e2", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
    set_double_default(lattice_element, "k2", 0.0);
    set_double_default(lattice_element, "h1", 0.0);
    set_double_default(lattice_element, "h2", 0.0);
    set_double_default(lattice_element, "hgap", 0.0);
    set_double_default(lattice_element, "fint", 0.0);
    set_double_default(lattice_element, "k3", 0.0);

}

Sbend_mad8_adaptor::~Sbend_mad8_adaptor()
{
}

Rbend_mad8_adaptor::Rbend_mad8_adaptor()
{
}

void
Rbend_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "angle", 0.0);
    set_double_default(lattice_element, "k1", 0.0);
    set_double_default(lattice_element, "e1", 0.0);
    set_double_default(lattice_element, "e2", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
    set_double_default(lattice_element, "k2", 0.0);
    set_double_default(lattice_element, "h1", 0.0);
    set_double_default(lattice_element, "h2", 0.0);
    set_double_default(lattice_element, "hgap", 0.0);
    set_double_default(lattice_element, "fint", 0.0);
    set_double_default(lattice_element, "k3", 0.0);
}

Rbend_mad8_adaptor::~Rbend_mad8_adaptor()
{
}

Quadrupole_mad8_adaptor::Quadrupole_mad8_adaptor()
{
}

void
Quadrupole_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "k1", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Quadrupole_mad8_adaptor::~Quadrupole_mad8_adaptor()
{
}

Sextupole_mad8_adaptor::Sextupole_mad8_adaptor()
{
}

void
Sextupole_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "k2", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Sextupole_mad8_adaptor::~Sextupole_mad8_adaptor()
{
}

Octupole_mad8_adaptor::Octupole_mad8_adaptor()
{
}

void
Octupole_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "k3", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Octupole_mad8_adaptor::~Octupole_mad8_adaptor()
{
}

Multipole_mad8_adaptor::Multipole_mad8_adaptor()
{
}

void
Multipole_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "lrad", 0.0);
    set_double_default(lattice_element, "k0l", 0.0);
    set_double_default(lattice_element, "t0", 0.0);
    set_double_default(lattice_element, "k1l", 0.0);
    set_double_default(lattice_element, "t1", 0.0);
    set_double_default(lattice_element, "k2l", 0.0);
    set_double_default(lattice_element, "t2", 0.0);
    set_double_default(lattice_element, "k3l", 0.0);
    set_double_default(lattice_element, "t3", 0.0);
    set_double_default(lattice_element, "k4l", 0.0);
    set_double_default(lattice_element, "t4", 0.0);
    set_double_default(lattice_element, "k5l", 0.0);
    set_double_default(lattice_element, "t5", 0.0);
    set_double_default(lattice_element, "k6l", 0.0);
    set_double_default(lattice_element, "t6", 0.0);
    set_double_default(lattice_element, "k7l", 0.0);
    set_double_default(lattice_element, "t7", 0.0);
    set_double_default(lattice_element, "k8l", 0.0);
    set_double_default(lattice_element, "t8", 0.0);
    set_double_default(lattice_element, "k9l", 0.0);
    set_double_default(lattice_element, "t9", 0.0);
}

Multipole_mad8_adaptor::~Multipole_mad8_adaptor()
{
}

Solenoid_mad8_adaptor::Solenoid_mad8_adaptor()
{
}

void
Solenoid_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "ks", 0.0);
}

Solenoid_mad8_adaptor::~Solenoid_mad8_adaptor()
{
}

Hkicker_mad8_adaptor::Hkicker_mad8_adaptor()
{
}

void
Hkicker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "kick", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Hkicker_mad8_adaptor::~Hkicker_mad8_adaptor()
{
}

Vkicker_mad8_adaptor::Vkicker_mad8_adaptor()
{
}

void
Vkicker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "kick", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Vkicker_mad8_adaptor::~Vkicker_mad8_adaptor()
{
}

Kicker_mad8_adaptor::Kicker_mad8_adaptor()
{
}

void
Kicker_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "hkick", 0.0);
    set_double_default(lattice_element, "vkick", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Kicker_mad8_adaptor::~Kicker_mad8_adaptor()
{
}

Rfcavity_mad8_adaptor::Rfcavity_mad8_adaptor()
{
}

void
Rfcavity_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "volt", 0.0);
    set_double_default(lattice_element, "lag", 0.0);
    set_double_default(lattice_element, "harmon", 0.0);
    set_double_default(lattice_element, "betrf", 0.0);
    set_double_default(lattice_element, "pg", 0.0);
    set_double_default(lattice_element, "shunt", 0.0);
    set_double_default(lattice_element, "tfill", 0.0);
}

Rfcavity_mad8_adaptor::~Rfcavity_mad8_adaptor()
{
}

Elseparator_mad8_adaptor::Elseparator_mad8_adaptor()
{
}

void
Elseparator_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "e", 0.0);
    set_double_default(lattice_element, "tilt", 0.0);
}

Elseparator_mad8_adaptor::~Elseparator_mad8_adaptor()
{
}

Hmonitor_mad8_adaptor::Hmonitor_mad8_adaptor()
{
}

void
Hmonitor_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
}

Hmonitor_mad8_adaptor::~Hmonitor_mad8_adaptor()
{
}

Vmonitor_mad8_adaptor::Vmonitor_mad8_adaptor()
{
}

void
Vmonitor_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
}

Vmonitor_mad8_adaptor::~Vmonitor_mad8_adaptor()
{
}

Monitor_mad8_adaptor::Monitor_mad8_adaptor()
{
}

void
Monitor_mad8_adaptor::set_default_attributes(Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
}

Monitor_mad8_adaptor::~Monitor_mad8_adaptor()
{
}

Instrument_mad8_adaptor::Instrument_mad8_adaptor()
{
}

void
Instrument_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
}

Instrument_mad8_adaptor::~Instrument_mad8_adaptor()
{
}

Ecollimator_mad8_adaptor::Ecollimator_mad8_adaptor()
{
}

void
Ecollimator_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "xsize", 0.0);
    set_double_default(lattice_element, "ysize", 0.0);
}

Ecollimator_mad8_adaptor::~Ecollimator_mad8_adaptor()
{
}

Rcollimator_mad8_adaptor::Rcollimator_mad8_adaptor()
{
}

void
Rcollimator_mad8_adaptor::set_default_attributes(
        Lattice_element & lattice_element)
{
    set_double_default(lattice_element, "l", 0.0);
    set_double_default(lattice_element, "xsize", 0.0);
    set_double_default(lattice_element, "ysize", 0.0);
}

Rcollimator_mad8_adaptor::~Rcollimator_mad8_adaptor()
{
}

