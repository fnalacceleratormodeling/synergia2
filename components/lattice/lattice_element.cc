#include "lattice_element.h"
#include <algorithm>
#include <iostream>

Lattice_element::Lattice_element(std::string const& type,
        std::string const& name) :
    type(type), name(name), ancestors(), double_attributes(),
            string_attributes(), length_attribute_name("l"),
            bend_angle_attribute_name("angle")
{

}

Lattice_element::Lattice_element(Lattice_element const& lattice_element) :
            type(lattice_element.type),
            name(lattice_element.name),
            ancestors(),
            double_attributes(),
            string_attributes(),
            length_attribute_name(lattice_element.length_attribute_name),
            bend_angle_attribute_name(lattice_element.bend_angle_attribute_name)
{
    std::copy(lattice_element.ancestors.begin(),
            lattice_element.ancestors.end(), std::inserter(ancestors,
                    ancestors.begin()));
    std::copy(lattice_element.double_attributes.begin(),
            lattice_element.double_attributes.end(), std::inserter(
                    double_attributes, double_attributes.begin()));
    std::copy(lattice_element.string_attributes.begin(),
            lattice_element.string_attributes.end(), std::inserter(
                    string_attributes, string_attributes.begin()));
}

std::string const &
Lattice_element::get_type() const
{
    return type;
}

std::string const &
Lattice_element::get_name() const
{
    return name;
}

void
Lattice_element::add_ancestor(std::string const& ancestor)
{
    ancestors.push_back(ancestor);
}

std::list<std::string > const&
Lattice_element::get_ancestors() const
{
    return ancestors;
}

void
Lattice_element::set_double_attribute(std::string const& name, double value)
{
    double_attributes[name] = value;
}

bool
Lattice_element::has_double_attribute(std::string const& name) const
{
    return (double_attributes.count(name) > 0);
}

double
Lattice_element::get_double_attribute(std::string const& name) const
{
    std::map<std::string, double >::const_iterator iter =
            double_attributes.find(name);
    return iter->second;
}

std::map<std::string, double > const &
Lattice_element::get_double_attributes() const
{
    return double_attributes;
}

void
Lattice_element::set_string_attribute(std::string const& name,
        std::string const& value)
{
    string_attributes[name] = value;
}

bool
Lattice_element::has_string_attribute(std::string const& name) const
{
    return (string_attributes.count(name) > 0);
}

std::string const&
Lattice_element::get_string_attribute(std::string const& name) const
{
    std::map<std::string, std::string >::const_iterator result =
            string_attributes.find(name);
    return result->second;
}

void
Lattice_element::set_length_attribute_name(std::string const& attribute_name)
{
    length_attribute_name = attribute_name;
}

void
Lattice_element::set_bend_angle_attribute_name(
        std::string const& attribute_name)
{
    bend_angle_attribute_name = attribute_name;
}

std::map<std::string, std::string > const &
Lattice_element::get_string_attributes() const
{
    return string_attributes;
}

double
Lattice_element::get_length() const
{
    std::map<std::string, double >::const_iterator iter =
            double_attributes.find(length_attribute_name);
    double retval = 0.0;
    if (iter != double_attributes.end()) {
        retval = iter->second;
    }
    return retval;
}

double
Lattice_element::get_bend_angle() const
{
    std::map<std::string, double >::const_iterator iter =
            double_attributes.find(bend_angle_attribute_name);
    double retval = 0.0;
    if (iter != double_attributes.end()) {
        retval = iter->second;
    }
    return retval;
}

void
Lattice_element::print() const
{
    for (std::list<std::string >::const_iterator it = ancestors.begin(); it
            != ancestors.end(); ++it) {
        std::cout << (*it) << ":";
    }
    std::cout << " ";
    std::cout << name << ": ";
    bool first_attr = true;
    for (std::map<std::string, double >::const_iterator it =
            double_attributes.begin(); it != double_attributes.end(); ++it) {
        if (first_attr) {
            first_attr = false;
        } else {
            std::cout << ", ";
        }
        std::cout << it->first << "=" << it->second;
    }
    for (std::map<std::string, std::string >::const_iterator it =
            string_attributes.begin(); it != string_attributes.end(); ++it) {
        if (first_attr) {
            first_attr = false;
        } else {
            std::cout << ", ";
        }
        std::cout << it->first << "=" << it->second;
    }
    std::cout << std::endl;
}
Set_default_attributes_fn_map
get_standard_default_attributes_fn_map()
{
    Set_default_attributes_fn_map map;
    map["marker"] = set_default_attributes_marker_mad8;
    map["drift"] = set_default_attributes_drift_mad8;
    map["sbend"] = set_default_attributes_sbend_mad8;
    map["rbend"] = set_default_attributes_rbend_mad8;
    map["quadrupole"] = set_default_attributes_quadrupole_mad8;
    map["sextupole"] = set_default_attributes_sextupole_mad8;
    map["octupole"] = set_default_attributes_octupole_mad8;
    map["multipole"] = set_default_attributes_multipole_mad8;
    map["solenoid"] = set_default_attributes_solenoid_mad8;
    map["hkicker"] = set_default_attributes_hkicker_mad8;
    map["vkicker"] = set_default_attributes_vkicker_mad8;
    map["kicker"] = set_default_attributes_kicker_mad8;
    map["rfcavity"] = set_default_attributes_rfcavity_mad8;
    map["elseparator"] = set_default_attributes_elseperator_mad8;
    map["hmonitor"] = set_default_attributes_hmonitor_mad8;
    map["vmonitor"] = set_default_attributes_vmonitor_mad8;
    map["monitor"] = set_default_attributes_monitor_mad8;
    map["instrument"] = set_default_attributes_instrument_mad8;
    map["ecollimator"] = set_default_attributes_ecollimator_mad8;
    map["rcollimator"] = set_default_attributes_rcollimator_mad8;

    return map;
}

void
set_default_attributes_marker_mad8(Lattice_element &lattice_element)
{
}

void
set_default_attributes_drift_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

void
set_default_attributes_sbend_mad8(Lattice_element &lattice_element)
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

void
set_default_attributes_rbend_mad8(Lattice_element &lattice_element)
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

void
set_default_attributes_quadrupole_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("k1", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

void
set_default_attributes_sextupole_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("k2", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

void
set_default_attributes_octupole_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("k3", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

void
set_default_attributes_multipole_mad8(Lattice_element &lattice_element)
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

void
set_default_attributes_solenoid_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("ks", 0.0);
}

void
set_default_attributes_hkicker_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("kick", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

void
set_default_attributes_vkicker_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("kick", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

void
set_default_attributes_kicker_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("hkick", 0.0);
    lattice_element.set_double_attribute("vkick", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

void
set_default_attributes_rfcavity_mad8(Lattice_element &lattice_element)
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

void
set_default_attributes_elseperator_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("e", 0.0);
    lattice_element.set_double_attribute("tilt", 0.0);
}

void
set_default_attributes_hmonitor_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

void
set_default_attributes_vmonitor_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

void
set_default_attributes_monitor_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

void
set_default_attributes_instrument_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
}

void
set_default_attributes_ecollimator_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("xsize", 0.0);
    lattice_element.set_double_attribute("ysize", 0.0);
}

void
set_default_attributes_rcollimator_mad8(Lattice_element &lattice_element)
{
    lattice_element.set_double_attribute("l", 0.0);
    lattice_element.set_double_attribute("xsize", 0.0);
    lattice_element.set_double_attribute("ysize", 0.0);
}

