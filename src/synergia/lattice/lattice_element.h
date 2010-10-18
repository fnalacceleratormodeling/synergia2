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

public:
    Lattice_element();
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
