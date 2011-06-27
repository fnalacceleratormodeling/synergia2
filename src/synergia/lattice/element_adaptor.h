#ifndef ELEMENT_ADAPTOR_H_
#define ELEMENT_ADAPTOR_H_

#include <map>
#include <string>
#include <list>

#include <boost/shared_ptr.hpp>
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/chef_elements.h"

#define THINPOLE 1

class Element_adaptor
{
public:
    Element_adaptor();
    void
    set_double_default(Lattice_element & lattice_element,
            std::string const& name, double value);
    void
    set_string_default(Lattice_element & lattice_element,
            std::string const& name, std::string const& value);
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    virtual Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Element_adaptor();
};

typedef boost::shared_ptr<Element_adaptor > Element_adaptor_sptr;

class Element_adaptor_map
{
private:
    std::map<std::string, Element_adaptor_sptr > adaptor_map;

public:
    Element_adaptor_map();
    void
    set_adaptor(std::string const& name,
            Element_adaptor_sptr element_adaptor_sptr);
    bool
    has_adaptor(std::string const& name) const;
    Element_adaptor_sptr
    get_adaptor(std::string const& name) const;
    std::list<std::string >
    get_adaptor_names() const;
    ~Element_adaptor_map();
};

typedef boost::shared_ptr<Element_adaptor_map > Element_adaptor_map_sptr;

class Marker_mad8_adaptor : public Element_adaptor
{
public:
    Marker_mad8_adaptor();
    void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Marker_mad8_adaptor();
};

class Drift_mad8_adaptor : public Element_adaptor
{
public:
    Drift_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Drift_mad8_adaptor();
};

class Sbend_mad8_adaptor : public Element_adaptor
{
public:
    Sbend_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Sbend_mad8_adaptor();
};

class Rbend_mad8_adaptor : public Element_adaptor
{
public:
    Rbend_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Rbend_mad8_adaptor();
};

class Quadrupole_mad8_adaptor : public Element_adaptor
{
public:
    Quadrupole_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Quadrupole_mad8_adaptor();
};

class Sextupole_mad8_adaptor : public Element_adaptor
{
public:
    Sextupole_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Sextupole_mad8_adaptor();
};

class Octupole_mad8_adaptor : public Element_adaptor
{
public:
    Octupole_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Octupole_mad8_adaptor();
};

class Multipole_mad8_adaptor : public Element_adaptor
{
public:
    Multipole_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Multipole_mad8_adaptor();
};

#ifdef THINPOLE
// thinpoles are an CHEF addon not found in MAD8
class Thinpole_mad8_adaptor : public Element_adaptor
{
 public:
  Thinpole_mad8_adaptor();
  virtual void
    set_default_attributes(Lattice_element & lattice_element);
  Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
  virtual
    ~Thinpole_mad8_adaptor();
};
#endif /* THINPOLE */

class Solenoid_mad8_adaptor : public Element_adaptor
{
public:
    Solenoid_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    virtual
    ~Solenoid_mad8_adaptor();
};

class Hkicker_mad8_adaptor : public Element_adaptor
{
public:
    Hkicker_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Hkicker_mad8_adaptor();
};

class Vkicker_mad8_adaptor : public Element_adaptor
{
public:
    Vkicker_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Vkicker_mad8_adaptor();
};

class Kicker_mad8_adaptor : public Element_adaptor
{
public:
    Kicker_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Kicker_mad8_adaptor();
};

class Rfcavity_mad8_adaptor : public Element_adaptor
{
public:
    Rfcavity_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Rfcavity_mad8_adaptor();
};

class Elseparator_mad8_adaptor : public Element_adaptor
{
public:
    Elseparator_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    virtual
    ~Elseparator_mad8_adaptor();
};

class Hmonitor_mad8_adaptor : public Element_adaptor
{
public:
    Hmonitor_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Hmonitor_mad8_adaptor();
};

class Vmonitor_mad8_adaptor : public Element_adaptor
{
public:
    Vmonitor_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Vmonitor_mad8_adaptor();
};

class Monitor_mad8_adaptor : public Element_adaptor
{
public:
    Monitor_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual
    ~Monitor_mad8_adaptor();
};

class Instrument_mad8_adaptor : public Element_adaptor
{
public:
    Instrument_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    virtual
    ~Instrument_mad8_adaptor();
};

class Ecollimator_mad8_adaptor : public Element_adaptor
{
public:
    Ecollimator_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    virtual
    ~Ecollimator_mad8_adaptor();
};

class Rcollimator_mad8_adaptor : public Element_adaptor
{
public:
    Rcollimator_mad8_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element);
    virtual
    ~Rcollimator_mad8_adaptor();
};


#endif /* ELEMENT_ADAPTOR_H_ */
