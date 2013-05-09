#ifndef ELEMENT_ADAPTOR_H_
#define ELEMENT_ADAPTOR_H_

#include <map>
#include <string>
#include <list>

#include <boost/shared_ptr.hpp>
#include "synergia/utils/serialization.h"

#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/chef_elements.h"

class Element_adaptor
{
    Lattice_element_sptr default_element_sptr;
public:
    Element_adaptor();
    Lattice_element_sptr
    get_default_element_sptr();
    Lattice_element &
    get_default_element();
    void
    set_double_default(Lattice_element & lattice_element,
            std::string const& name, double value);
    void
    set_string_default(Lattice_element & lattice_element,
            std::string const& name, std::string const& value);
    virtual void
    set_defaults(Lattice_element & lattice_element);
    virtual void
    set_derived_attributes_internal(Lattice_element & lattice_element);
    virtual void
    set_derived_attributes_external(Lattice_element & lattice_element,
            double lattice_length, double beta);
    virtual Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Element_adaptor();
};

typedef boost::shared_ptr<Element_adaptor > Element_adaptor_sptr; // syndoc:include

#endif /* ELEMENT_ADAPTOR_H_ */
