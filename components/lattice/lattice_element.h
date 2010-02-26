#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <string>
#include <map>
#include <list>

class Lattice_element;
typedef void
(*set_default_attributes_fn)(Lattice_element &);
typedef std::map<std::string, set_default_attributes_fn >
        Set_default_attributes_fn_map;

class Lattice_element
{
private:
    std::string type;
    std::string name;
    std::list<std::string > ancestors;
    std::map<std::string, double > double_attributes;
    std::map<std::string, std::string > string_attributes;
    std::string length_attribute_name;
    std::string bend_angle_attribute_name;

    void
    set_default_attributes(Set_default_attributes_fn_map const& map);

public:
    Lattice_element(std::string const& type, std::string const& name,
            Set_default_attributes_fn_map const& map);
    Lattice_element(std::string const& type, std::string const& name);
    Lattice_element(Lattice_element const& lattice_element);
    std::string const &
    get_type() const;
    std::string const &
    get_name() const;
    void
    add_ancestor(std::string const& ancestor);
    std::list<std::string > const&
    get_ancestors() const;
    void
    set_double_attribute(std::string const& name, double value);
    bool
    has_double_attribute(std::string const& name) const;
    double
    get_double_attribute(std::string const& name) const;
    std::map<std::string, double > const &
    get_double_attributes() const;
    void
    set_string_attribute(std::string const& name, std::string const& value);
    bool
    has_string_attribute(std::string const& name) const;
    std::string const&
    get_string_attribute(std::string const& name) const;
    void
    set_length_attribute_name(std::string const& attribute_name);
    void
    set_bend_angle_attribute_name(std::string const& attribute_name);
    std::map<std::string, std::string > const &
    get_string_attributes() const;
    double
    get_length() const;
    double
    get_bend_angle() const;
};

Set_default_attributes_fn_map
get_standard_default_attributes_fn_map();

void
set_default_attributes_marker_mad8(Lattice_element &lattice_element);
void
set_default_attributes_drift_mad8(Lattice_element &lattice_element);
void
set_default_attributes_sbend_mad8(Lattice_element &lattice_element);
void
set_default_attributes_rbend_mad8(Lattice_element &lattice_element);
void
set_default_attributes_quadrupole_mad8(Lattice_element &lattice_element);
void
set_default_attributes_sextupole_mad8(Lattice_element &lattice_element);
void
set_default_attributes_octupole_mad8(Lattice_element &lattice_element);
void
set_default_attributes_multipole_mad8(Lattice_element &lattice_element);
void
set_default_attributes_solenoid_mad8(Lattice_element &lattice_element);
void
set_default_attributes_hkicker_mad8(Lattice_element &lattice_element);
void
set_default_attributes_vkicker_mad8(Lattice_element &lattice_element);
void
set_default_attributes_kicker_mad8(Lattice_element &lattice_element);
void
set_default_attributes_rfcavity_mad8(Lattice_element &lattice_element);
void
set_default_attributes_elseperator_mad8(Lattice_element &lattice_element);
void
set_default_attributes_hmonitor_mad8(Lattice_element &lattice_element);
void
set_default_attributes_vmonitor_mad8(Lattice_element &lattice_element);
void
set_default_attributes_monitor_mad8(Lattice_element &lattice_element);
void
set_default_attributes_instrument_mad8(Lattice_element &lattice_element);
void
set_default_attributes_ecollimator_mad8(Lattice_element &lattice_element);
void
set_default_attributes_rcollimator_mad8(Lattice_element &lattice_element);

#endif /* LATTICE_ELEMENT_H_ */
