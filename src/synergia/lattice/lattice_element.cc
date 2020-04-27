#include "lattice_element.h"
#include "lattice.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <sstream>


namespace
{
    using type_map_t = std::map<std::string, element_type>;

    const type_map_t type_map = { 
        { element_type_name::generic,    element_type::generic},
        { element_type_name::drift,      element_type::drift},
        { element_type_name::rbend,      element_type::rbend},
        { element_type_name::sbend,      element_type::sbend},
        { element_type_name::quadrupole, element_type::quadrupole},
        { element_type_name::multipole,  element_type::multipole},
        { element_type_name::rfcavity,   element_type::rfcavity},
        { element_type_name::hkicker,    element_type::hkicker},
        { element_type_name::vkicker,    element_type::vkicker},
        { element_type_name::kicker,     element_type::kicker},
        { element_type_name::monitor,    element_type::monitor},
        { element_type_name::sextupole,  element_type::sextupole},
        { element_type_name::octupole,   element_type::octupole},
        { element_type_name::marker,     element_type::marker},
        { element_type_name::instrument, element_type::instrument},
        { element_type_name::rcollimator,element_type::rcollimator},
        { element_type_name::nllens,     element_type::nllens},
    };

    element_type find_type(std::string const & stype)
    {
        auto r = type_map.find(stype);
        if (r == type_map.end()) throw std::runtime_error("invalid element type " + stype);
        return r->second;
    }
}


Lattice_element::Lattice_element() 
    : name("")
    , stype("generic")
    , type(find_type(stype))
    , ancestors()
    , double_attributes()
    , string_attributes()
    , vector_attributes()
    , length_attribute_name("l")
    , bend_angle_attribute_name("angle")
    , revision(0)
    , lattice_ptr(nullptr)
{
}

Lattice_element::Lattice_element(
        std::string const & type,
        std::string const & name) 
    : name(name)
    , stype(type)
    , type(find_type(stype))
    , ancestors()
    , double_attributes()
    , string_attributes()
    , vector_attributes()
    , length_attribute_name("l")
    , bend_angle_attribute_name("angle")
    , revision(0)
    , lattice_ptr(nullptr)
{
}

Lattice_element::Lattice_element(Lsexpr const & lsexpr)
    : name("")
    , stype("generic")
    , type(find_type(stype))
    , ancestors()
    , double_attributes()
    , string_attributes()
    , vector_attributes()
    , length_attribute_name("l")
    , bend_angle_attribute_name("angle")
    , revision(0)
    , lattice_ptr(nullptr)
{
    for (auto const& lse : lsexpr)
    {
        if (lse.is_labeled()) 
        {
            if (lse.get_label() == "type") 
            {
                stype = lse.get_string();
                type = find_type(stype);
            } 
            else if (lse.get_label() == "name") 
            {
                name = lse.get_string();
            } 
            else if (lse.get_label() == "ancestors") 
            {
                auto ancestors_vector = lse.get_string_vector();
                std::copy(ancestors_vector.begin(), ancestors_vector.end(),
                          std::back_inserter(ancestors));
            } 
            else if (lse.get_label() == "double_attributes") 
            {
                for (auto const& attr : lse)
                    double_attributes[attr.get_label()] = attr.get_double();
            } 
            else if (lse.get_label() == "string_attributes") 
            {
                for (auto const& attr : lse)
                    string_attributes[attr.get_label()] = attr.get_string();
            } 
            else if (lse.get_label() == "vector_attributes") 
            {
                for (auto const& attr : lse)
                    vector_attributes[attr.get_label()] = attr.get_double_vector();
            }
        } 
        else 
        {
            if (!lse.is_atomic())
            {
                for (auto const& attr : lse)
                    double_attributes[attr.get_label()] = attr.get_double();
            }
        }
    }
}

Lsexpr
Lattice_element::as_lsexpr() const
{
    Lsexpr retval;
#if 0
    retval.push_back(Lsexpr(type, "type"));
    retval.push_back(Lsexpr(name, "name"));
    if (double_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, double>::const_iterator it =
                 double_attributes.begin();
             it != double_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        if (!((string_attributes.size() == 0) &&
              (vector_attributes.size() == 0))) {
            attrs.set_label("double_attributes");
        }
    }
    if (double_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, double>::const_iterator it =
                 double_attributes.begin();
             it != double_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        if (!((string_attributes.size() == 0) &&
              (vector_attributes.size() == 0))) {
            attrs.set_label("double_attributes");
        }
        retval.push_back(attrs);
    }
    if (string_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, std::string>::const_iterator it =
                 string_attributes.begin();
             it != string_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        attrs.set_label("string_attributes");
        retval.push_back(attrs);
    }
    if (vector_attributes.size() > 0) {
        Lsexpr attrs;
        for (std::map<std::string, std::vector<double> >::const_iterator it =
                 vector_attributes.begin();
             it != vector_attributes.end(); ++it) {
            Lsexpr attr(it->second);
            attr.set_label(it->first);
            attrs.push_back(attr);
        }
        attrs.set_label("vector_attributes");
        retval.push_back(attrs);
    }
    if (ancestors.size() > 0) {
        std::vector<std::string> ancestors_vector;
        std::copy(ancestors.begin(), ancestors.end(),
                  std::back_inserter(ancestors_vector));
        retval.push_back(Lsexpr(ancestors_vector, "ancestors"));
    }
#endif
    return retval;
}

std::string const &
Lattice_element::get_type_name() const
{
    return stype;
}

element_type
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
Lattice_element::set_double_attribute(
        std::string const & name, 
        double value,
        bool increment_revision)
{
    double_attributes[name] = value;
    if (increment_revision) ++revision;
}

void
Lattice_element::set_default_double_attribute(
        std::string const & name, 
        double value,
        bool increment_revision )
{
    if (!has_double_attribute(name))
    {
        double_attributes[name] = value;
        if (increment_revision) ++revision;
    }
}

bool
Lattice_element::has_double_attribute(std::string const & name) const
{
    bool retval = (double_attributes.count(name) > 0);
    return retval;
}

double
Lattice_element::get_double_attribute(std::string const& name) const
{
    auto r = double_attributes.find(name);
    if (r == double_attributes.end()) 
    { 
        throw std::runtime_error( 
                "Lattice_element::get_double_attribute: element "
                + this->name + " of type " + stype
                + " has no double attribute '" + name + "'");
    } 
    else 
    {
        return r->second;
    }
}

double
Lattice_element::get_double_attribute(std::string const& name, double val) const
{
    auto r = double_attributes.find(name);
    if (r == double_attributes.end()) 
    {
        return val;
    } 
    else 
    {
        return r->second;
    }
}

void
Lattice_element::set_string_attribute(
        std::string const & name,
        std::string const & value, 
        bool increment_revision)
{
    string_attributes[name] = value;
    if (increment_revision) ++revision;
}

void
Lattice_element::set_default_string_attribute(
        std::string const & name,
        std::string const & value, 
        bool increment_revision)
{
    if (!has_string_attribute(name))
    {
        string_attributes[name] = value;
        if (increment_revision) ++revision;
    }
}

bool
Lattice_element::has_string_attribute(std::string const & name) const
{
    bool retval = (string_attributes.count(name) > 0);
    return retval;
}

std::string const&
Lattice_element::get_string_attribute(std::string const & name) const
{
    auto r = string_attributes.find(name);
    if (r == string_attributes.end()) 
    { 
        throw std::runtime_error( 
                "Lattice_element::get_string_attribute: element "
                + this->name + " of type " + stype
                + " has no string attribute '" + name + "'");
    } 
    else 
    {
        return r->second;
    }
}

std::string const&
Lattice_element::get_string_attribute(std::string const & name, std::string const & val) const
{
    auto r = string_attributes.find(name);
    if (r == string_attributes.end()) return val;
    else return r->second;
}

void
Lattice_element::set_vector_attribute(
        std::string const & name,
        std::vector<double> const & value, 
        bool increment_revision)
{
    vector_attributes[name] = value;
    if (increment_revision) ++revision;
}

bool
Lattice_element::has_vector_attribute(std::string const & name) const
{
    bool retval = (vector_attributes.count(name) > 0);
    return retval;
}

std::vector<double> const &
Lattice_element::get_vector_attribute(std::string const & name) const
{
    auto r = vector_attributes.find(name);
    if (r == vector_attributes.end()) 
    { 
        throw std::runtime_error( 
                "Lattice_element::get_vector_attribute: element "
                + this->name + " of type " + stype
                + " has no vector attribute '" + name + "'");
    } 
    else 
    {
        return r->second;
    }
}

std::vector<double> const&
Lattice_element::get_vector_attribute(
        std::string const & name, 
        std::vector<double> const & val) const
{
    auto r = vector_attributes.find(name);
    if (r == vector_attributes.end()) return val;
    else return r->second;
}

void
Lattice_element::set_length_attribute_name(std::string const & attribute_name)
{
    length_attribute_name = attribute_name;
}

void
Lattice_element::set_bend_angle_attribute_name(std::string const & attribute_name)
{
    bend_angle_attribute_name = attribute_name;
}

double
Lattice_element::get_length() const
{
    return get_double_attribute(length_attribute_name, 0.0);
}

double
Lattice_element::get_bend_angle() const
{
    return get_double_attribute(bend_angle_attribute_name, 0.0);
}

long int
Lattice_element::get_revision() const
{
    return revision;
}

bool
Lattice_element::has_lattice() const
{
    return (lattice_ptr != 0);
}

void
Lattice_element::set_lattice(Lattice & lattice)
{
    lattice_ptr = &lattice;
}

Lattice const&
Lattice_element::get_lattice() const
{
    if(! has_lattice()) {
        throw std::runtime_error(
                    "Lattice_element::get_lattice: element not part of any lattice");
    }
    return *lattice_ptr;
}


std::string
Lattice_element::as_string() const
{
    std::stringstream sstream;
    sstream.precision(15);
    for (std::list<std::string >::const_iterator it = ancestors.begin();
            it != ancestors.end(); ++it) {
        sstream << (*it) << ":";
    }
    sstream << " " << stype << " ";
    sstream << name << ": ";
    bool first_attr = true;
    for (std::map<std::string, double >::const_iterator it =
            double_attributes.begin(); it != double_attributes.end(); ++it) {
        if (first_attr) {
            first_attr = false;
        } else {
            sstream << ", ";
        }
        sstream << it->first << "=" << it->second;
    }
    for (std::map<std::string, std::string >::const_iterator it =
            string_attributes.begin(); it != string_attributes.end(); ++it) {
        if (first_attr) {
            first_attr = false;
        } else {
            sstream << ", ";
        }
        sstream << it->first << "=" << it->second;
    }
    for (std::map<std::string, std::vector<double > >::const_iterator it =
             vector_attributes.begin(); it != vector_attributes.end(); ++it) {
        if (first_attr) {
            first_attr = false;
        } else {
            sstream << ", ";
        }
        sstream << it->first << "=" << "{";
        for (size_t i=0; i != (it->second).size(); ++i) {
            if (i) {
                sstream << ", ";
            }
            sstream << (it->second)[i];
        }
        sstream << "}";
    }
    return sstream.str();
}


void
Lattice_element::print() const
{
    std::cout << as_string() << std::endl;
}

