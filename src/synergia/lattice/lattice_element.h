#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <string>
#include <map>
#include <list>

#include <boost/shared_ptr.hpp>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/version.hpp>

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
    std::list<std::string > ancestors;
    std::map<std::string, double > double_attributes;
    std::map<std::string, std::string > string_attributes;
    std::string length_attribute_name;
    std::string bend_angle_attribute_name;
    long int revision;

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
    void
    set_double_attribute(std::string const& name, double value);

    /// Check for the existence of the named double attribute
    /// @param name attribute name
    bool
    has_double_attribute(std::string const& name) const;

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
    void
    set_string_attribute(std::string const& name, std::string const& value);

    /// Check for the existence of the named string attribute
    /// @param name attribute name
    bool
    has_string_attribute(std::string const& name) const;

    /// Get the value of the named string attribute
    /// @param name attribute name
    std::string const&
    get_string_attribute(std::string const& name) const;

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

    /// Get the entire dictionary of string attributes
    std::map<std::string, std::string > const &
    get_string_attributes() const;

    /// Get the Lattice_element's length
    double
    get_length() const;

    /// Get the Lattice_element's bend angle
    double
    get_bend_angle() const;

    /// Get the Lattice_element's revision number
    long int
    get_revision() const;

    /// Increment the Lattice_element's revision number
    void
    increment_revision();

    /// Print a human-readable description of the Lattice_element
    /// The Python version of the function is named "print_".
    void
    print() const;

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(type) & BOOST_SERIALIZATION_NVP(name)
                    & BOOST_SERIALIZATION_NVP(ancestors)
                    & BOOST_SERIALIZATION_NVP(double_attributes)
                    & BOOST_SERIALIZATION_NVP(string_attributes)
                    & BOOST_SERIALIZATION_NVP(length_attribute_name)
                    & BOOST_SERIALIZATION_NVP(bend_angle_attribute_name);
        }
};

typedef boost::shared_ptr<Lattice_element > Lattice_element_sptr;
typedef std::list<Lattice_element_sptr > Lattice_elements;

#endif /* LATTICE_ELEMENT_H_ */
