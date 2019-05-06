#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <string>
#include <map>
#include <list>
#include <vector>

#include "synergia/utils/cereal.h"
#include "synergia/utils/lsexpr.h"


enum class element_type
{
    generic,

    drift,
    rbend,
    sbend,
    quadrupole,
};

namespace element_type_name
{
    constexpr char const * generic = "generic";
    constexpr char const * drift   = "drift";
}

class Lattice_element;
typedef std::shared_ptr<Lattice_element> Lattice_element_sptr;

class Lattice_data;

/// The Lattice_element class contains the description of a single
/// lattice element. Each element has a name, a (string) type and
/// dictionaries of named double and string attributes.
/// Lattice structure is described by a list of ancestors stored in
/// an element.
class Lattice_element
{

private:

    std::string  name;

    std::string  stype;
    element_type type;

    std::list<std::string > ancestors;

    std::map<std::string, double> double_attributes;
    std::map<std::string, std::string> string_attributes;
    std::map<std::string, std::vector<double>> vector_attributes;

    std::string length_attribute_name;
    std::string bend_angle_attribute_name;

    long int revision;

    Lattice_data * lattice_ptr;

public:

    /// Construct a Lattice_element with an empty name and type.
    Lattice_element();

    /// Construct a Lattice_element.
    /// @param name name
    /// @param type type
    Lattice_element(std::string const & type, std::string const & name);

    /// Construct a Lattice_element from the Lsexpr representation
    /// @param lsexpr representation
    explicit Lattice_element(Lsexpr const& lsexpr);

    /// Extract an Lsexpr representation of the Lattice_element
    Lsexpr
    as_lsexpr() const;

    /// Get the type
    std::string const &
    get_type_name() const;

    element_type
    get_type() const;

    /// Get the name
    std::string const &
    get_name() const;

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
    set_double_attribute(
            std::string const & name, 
            double value,
            bool increment_revision = true);

    void
    set_default_double_attribute(
            std::string const & name, 
            double value,
            bool increment_revision = true);

    /// Check for the existence of the named double attribute
    /// @param name attribute name
    bool
    has_double_attribute(std::string const & name) const;

    /// Get the value of the named double attribute
    /// @param name attribute name
    double
    get_double_attribute(std::string const & name) const;

    /// Get the value of the named double attribute
    /// @param name attribute name
    /// @param val default value if the specified attribute doesnt exist
    double
    get_double_attribute(std::string const & name, double val) const;

    /// Set the value of the named string attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_string_attribute(
            std::string const & name, 
            std::string const & value,
            bool increment_revision = true);

    void
    set_default_string_attribute(
            std::string const & name, 
            std::string const & value,
            bool incremnt_revision = true);

    /// Check for the existence of the named string attribute
    /// @param name attribute name
    bool
    has_string_attribute(std::string const & name) const;

    /// Get the value of the named string attribute
    /// @param name attribute name
    std::string const&
    get_string_attribute(std::string const & name) const;

    /// Get the value of the named string attribute
    /// @param name attribute name
    /// @param val default value if the specified attribute doesnt exist
    std::string const&
    get_string_attribute(std::string const & name, std::string const & val) const;

    /// Set the value of the named vector attribute
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_vector_attribute(
            std::string const & name,
            std::vector<double> const & value, 
            bool increment_revision = true);

    /// Check for the existence of the named vector attribute
    /// @param name attribute name
    bool
    has_vector_attribute(std::string const & name) const;

    /// Get the value of the named vector attribute
    /// @param name attribute name
    std::vector<double> const &
    get_vector_attribute(std::string const & name) const;

    /// Get the value of the named vector attribute
    /// @param name attribute name
    /// @param val default value if the specified attribute doesnt exist
    std::vector<double> const &
    get_vector_attribute(std::string const & name, std::vector<double> const & val) const;

    /// Set the attribute name to be used to determine the length
    /// of the Lattice_element
    /// @param attribute_name attribute name
    void
    set_length_attribute_name(std::string const & attribute_name);

    /// Set the attribute name to be used to determine the bend_angle
    /// of the Lattice_element
    /// @param attribute_name attribute name
    void
    set_bend_angle_attribute_name(std::string const & attribute_name);

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
    set_lattice(Lattice_data & lattice);

    /// Get a reference to the parent lattice
    Lattice_data &
    get_lattice();

    /// Get a reference to the parent lattice
    Lattice_data const &
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

//typedef std::list<Lattice_element_sptr > Lattice_elements; // syndoc:include
typedef std::list<Lattice_element> Lattice_elements; // syndoc:include

#endif /* LATTICE_ELEMENT_H_ */
