#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <string>
#include <map>
#include <list>
#include <vector>

#include <boost/shared_ptr.hpp>
#include "synergia/utils/serialization.h"

class Lattice_element;
typedef boost::shared_ptr<Lattice_element > Lattice_element_sptr; // syndoc:include

class Lattice;

/// The Lattice_element class contains the description of a single
/// lattice element. Each element has a name, a (string) type and
/// dictionaries of named double and string attributes.
/// Lattice structure is described by a list of ancestors stored in
/// an element.
class Lattice_element
{
private:
    std::string type;
    std::string name;
    Lattice_element_sptr default_element_sptr;
    std::list<std::string > ancestors;
    std::map<std::string, double > double_attributes;
    std::map<std::string, std::string > string_attributes;
    std::map<std::string, std::vector<double > > vector_attributes;
    std::string length_attribute_name;
    std::string bend_angle_attribute_name;
    long int revision;
    bool needs_internal_derive, needs_external_derive;
    Lattice *lattice_ptr;

public:
    /// Construct a Lattice_element with an empty name and type.
    Lattice_element();

    /// Construct a Lattice_element.
    /// @param name name
    /// @param type type
    Lattice_element(std::string const& type, std::string const& name);

    /// Copy constructor.
    Lattice_element(Lattice_element const& lattice_element);

    /// Get the type
    std::string const &
    get_type() const;

    /// Get the name
    std::string const &
    get_name() const;

    /// Set the defaults for this element
    void
    set_default_element(Lattice_element_sptr default_element_sptr);

    /// Add an ancestor to the list of ancestors
    /// @param ancestor ancestor name
    void
    add_ancestor(std::string const& ancestor);

    /// Get the list of ancestors
    std::list<std::string > const&
    get_ancestors() const;

    /// Set the value of the named double attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_double_attribute(std::string const& name, double value,
            bool increment_revision = true);

    /// Check for the existence of the named double attribute
    /// @param name attribute name
    bool
    has_double_attribute(std::string const& name,
            bool include_default = true) const;

    /// Get the value of the named double attribute
    /// @param name attribute name
    double
    get_double_attribute(std::string const& name) const;

    /// Get the entire dictionary of double attributes
    std::map<std::string, double > const &
    get_double_attributes() const;

    /// Set the value of the named string attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_string_attribute(std::string const& name, std::string const& value,
            bool increment_revision = true);

    /// Check for the existence of the named string attribute
    /// @param name attribute name
    bool
    has_string_attribute(std::string const& name,
            bool include_default = true) const;

    /// Get the value of the named string attribute
    /// @param name attribute name
    std::string const&
    get_string_attribute(std::string const& name) const;

    /// Get the entire dictionary of string attributes
    std::map<std::string, std::string > const &
    get_string_attributes() const;

    /// Set the value of the named vector attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_vector_attribute(std::string const& name,
            std::vector<double > const& value, bool increment_revision = true);

    /// Check for the existence of the named vector attribute
    /// @param name attribute name
    bool
    has_vector_attribute(std::string const& name,
            bool include_default = true) const;

    /// Get the value of the named vector attribute
    /// @param name attribute name
    std::vector<double > const&
    get_vector_attribute(std::string const& name) const;

    /// Get the entire dictionary of vector attributes
    std::map<std::string, std::vector<double > > const &
    get_vector_attributes() const;

    /// Set the attribute name to be used to determine the length
    /// of the Lattice_element
    /// @param attribute_name attribute name
    void
    set_length_attribute_name(std::string const& attribute_name);

    /// Set the attribute name to be used to determine the bend_angle
    /// of the Lattice_element
    /// @param attribute_name attribute name
    void
    set_bend_angle_attribute_name(std::string const& attribute_name);

    /// Set whether the element needs to determine some of its parameters
    /// from its other parameters
    void
    set_needs_internal_derive(bool value);

    /// Get whether the element needs to determine some of its parameters
    /// from its other parameters
    bool
    get_needs_internal_derive() const;

    /// Set whether the element needs to determine some of its parameters
    /// from the lattice length and/or reference particle
    void
    set_needs_external_derive(bool value);

    /// Get whether the element needs to determine some of its parameters
    /// from the lattice length and/or reference particle
    bool
    get_needs_external_derive() const;

    /// Get the Lattice_element's length
    double
    get_length() const;

    /// Get the Lattice_element's bend angle
    double
    get_bend_angle() const;

    /// Get the Lattice_element's revision number
    long int
    get_revision() const;

    /// Check whether the element has a reference to a parent lattice
    bool
    has_lattice() const;

    /// Set the reference to the parent lattice
    void
    set_lattice(Lattice &lattice);

    /// Get a reference to the parent lattice
    Lattice &
    get_lattice();

    /// Get a reference to the parent lattice
    Lattice &
    get_lattice() const;

    /// Return a human-readable description of the Lattice_element
    std::string
    as_string() const;

    /// Print a human-readable description of the Lattice_element
    /// The Python version of the function is named "print_".
    void
    print() const;
    
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

typedef std::list<Lattice_element_sptr > Lattice_elements; // syndoc:include

#endif /* LATTICE_ELEMENT_H_ */
