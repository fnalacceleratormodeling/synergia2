#include "lattice_element.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <sstream>

Lattice_element::Lattice_element() :
        type(""), name(""), default_element_sptr(), ancestors(), double_attributes(), string_attributes(), length_attribute_name(
                "l"), bend_angle_attribute_name("angle"), revision(0), needs_internal_derive(
                false), needs_external_derive(false)
{

}

Lattice_element::Lattice_element(std::string const& type,
        std::string const& name) :
        type(type), name(name), default_element_sptr(), ancestors(), double_attributes(), string_attributes(), length_attribute_name(
                "l"), bend_angle_attribute_name("angle"), revision(0), needs_internal_derive(
                false), needs_external_derive(false)
{

}

Lattice_element::Lattice_element(Lattice_element const& lattice_element) :
        type(lattice_element.type), name(lattice_element.name), default_element_sptr(
                lattice_element.default_element_sptr), ancestors(), double_attributes(), string_attributes(), length_attribute_name(
                lattice_element.length_attribute_name), bend_angle_attribute_name(
                lattice_element.bend_angle_attribute_name), revision(0), needs_internal_derive(
                lattice_element.needs_internal_derive), needs_external_derive(
                lattice_element.needs_external_derive)
{
    std::copy(lattice_element.ancestors.begin(),
            lattice_element.ancestors.end(),
            std::inserter(ancestors, ancestors.begin()));
    std::copy(lattice_element.double_attributes.begin(),
            lattice_element.double_attributes.end(),
            std::inserter(double_attributes, double_attributes.begin()));
    std::copy(lattice_element.string_attributes.begin(),
            lattice_element.string_attributes.end(),
            std::inserter(string_attributes, string_attributes.begin()));
    std::copy(lattice_element.vector_attributes.begin(),
            lattice_element.vector_attributes.end(),
            std::inserter(vector_attributes, vector_attributes.begin()));
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
Lattice_element::set_default_element(Lattice_element_sptr default_element_sptr)
{
    this->default_element_sptr = default_element_sptr;
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
Lattice_element::set_double_attribute(std::string const& name, double value,
        bool increment_revision)
{
    double_attributes[name] = value;
    if (increment_revision) {
        ++revision;
    }
}

bool
Lattice_element::has_double_attribute(std::string const& name,
        bool include_default) const
{
    bool retval = (double_attributes.count(name) > 0);
    if ((!retval) && include_default && default_element_sptr) {
        retval = default_element_sptr->has_double_attribute(name, false);
    }
    return retval;
}

double
Lattice_element::get_double_attribute(std::string const& name) const
{
    std::map<std::string, double >::const_iterator result =
            double_attributes.find(name);
    if (result == double_attributes.end()) {
        if (default_element_sptr
                && default_element_sptr->has_double_attribute(name, false)) {
            return default_element_sptr->get_double_attribute(name);
        } else {
            throw std::runtime_error(
                    "Lattice_element::get_double_attribute: element "
                            + this->name + " of type " + type
                            + " has no double attribute '" + name + "'");
        }
    } else {
        return result->second;
    }
}

std::map<std::string, double > const &
Lattice_element::get_double_attributes() const
{
    return double_attributes;
}

void
Lattice_element::set_string_attribute(std::string const& name,
        std::string const& value, bool increment_revision)
{
    string_attributes[name] = value;
    if (increment_revision) {
        ++revision;
    }
}

bool
Lattice_element::has_string_attribute(std::string const& name,
        bool include_default) const
{
    bool retval = (string_attributes.count(name) > 0);
    if ((!retval) && include_default && default_element_sptr) {
        retval = default_element_sptr->has_string_attribute(name, false);
    }
    return retval;
}

std::string const&
Lattice_element::get_string_attribute(std::string const& name) const
{
    std::map<std::string, std::string >::const_iterator result =
            string_attributes.find(name);
    if (result == string_attributes.end()) {
        if (default_element_sptr
                && default_element_sptr->has_string_attribute(name, false)) {
            return default_element_sptr->get_string_attribute(name);
        } else {
            throw std::runtime_error(
                    "Lattice_element::get_string_attribute: element "
                            + this->name + " of type " + type
                            + " has no string attribute '" + name + "'");
        }
    } else {
        return result->second;
    }
}

std::map<std::string, std::string > const &
Lattice_element::get_string_attributes() const
{
    return string_attributes;
}

void
Lattice_element::set_vector_attribute(std::string const& name,
        std::vector<double > const& value, bool increment_revision)
{
    vector_attributes[name] = value;
    if (increment_revision) {
        ++revision;
    }
}

bool
Lattice_element::has_vector_attribute(std::string const& name,
        bool include_default) const
{
    bool retval = (vector_attributes.count(name) > 0);
    if ((!retval) && include_default && default_element_sptr) {
        retval = default_element_sptr->has_vector_attribute(name, false);
    }
    return retval;
}

std::vector<double > const&
Lattice_element::get_vector_attribute(std::string const& name) const
{
    std::map<std::string, std::vector<double > >::const_iterator result =
            vector_attributes.find(name);
    if (result == vector_attributes.end()) {
        if (default_element_sptr
                && default_element_sptr->has_vector_attribute(name, false)) {
            return default_element_sptr->get_vector_attribute(name);
        } else {
            throw std::runtime_error(
                    "Lattice_element::get_vector_attribute: element "
                            + this->name + " of type " + type
                            + " has no vector attribute '" + name + "'");
        }
    } else {
        return result->second;
    }
}

std::map<std::string, std::vector<double > > const &
Lattice_element::get_vector_attributes() const
{
    return vector_attributes;
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

void
Lattice_element::set_needs_internal_derive(bool value)
{
    needs_internal_derive = value;
}

bool
Lattice_element::get_needs_internal_derive() const
{
    return needs_internal_derive;
}

void
Lattice_element::set_needs_external_derive(bool value)
{
    needs_external_derive = value;
}

bool
Lattice_element::get_needs_external_derive() const
{
    return needs_external_derive;
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

long int
Lattice_element::get_revision() const
{
    return revision;
}

std::string
Lattice_element::as_string() const
{
    std::stringstream sstream;
    for (std::list<std::string >::const_iterator it = ancestors.begin();
            it != ancestors.end(); ++it) {
        sstream << (*it) << ":";
    }
    sstream << " " << type << " ";
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


template<class Archive>
    void
    Lattice_element::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(type)
        & BOOST_SERIALIZATION_NVP(name)
        & BOOST_SERIALIZATION_NVP(default_element_sptr)
        & BOOST_SERIALIZATION_NVP(ancestors)
        & BOOST_SERIALIZATION_NVP(double_attributes)
        & BOOST_SERIALIZATION_NVP(string_attributes)
        & BOOST_SERIALIZATION_NVP(vector_attributes)
        & BOOST_SERIALIZATION_NVP(length_attribute_name)
        & BOOST_SERIALIZATION_NVP(bend_angle_attribute_name)
        & BOOST_SERIALIZATION_NVP(revision)
        & BOOST_SERIALIZATION_NVP(needs_internal_derive)
        & BOOST_SERIALIZATION_NVP(needs_external_derive);
    }

template
void
Lattice_element::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lattice_element::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Lattice_element::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lattice_element::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
