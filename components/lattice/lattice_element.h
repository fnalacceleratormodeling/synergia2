#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <string>
#include <map>

class Lattice_element
{
private:
    std::map<std::string, double > double_attributes;
    std::map<std::string, std::string > string_attributes;
    std::string type;
    std::string name;

public:
    Lattice_element(std::string const& type, std::string const& name);
    std::string const &
    get_type() const;
    std::string const &
    get_name() const;
    void
    set_double_attribute(std::string const& name, double value);
    bool
    has_double_attribute(std::string const& name) const;
    double
    get_double_attribute(std::string const& name);
    //	double length();
    //	double quadrupole_strength();
};

#endif /* LATTICE_ELEMENT_H_ */
