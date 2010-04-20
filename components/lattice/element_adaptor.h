#ifndef ELEMENT_ADAPTOR_H_
#define ELEMENT_ADAPTOR_H_

#include <map>
#include <string>
#include <list>

#include <boost/shared_ptr.hpp>
#include "components/lattice/lattice_element.h"

class Element_adaptor
{
public:
    Element_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element) = 0;
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
            Element_adaptor_sptr const& element_adaptor_sptr);
    Element_adaptor_sptr &
    get_adaptor(std::string const& name);
    std::list<std::string >
    get_adaptor_names() const;
    ~Element_adaptor_map();
};

typedef boost::shared_ptr<Element_adaptor_map > Element_adaptor_map_sptr;

class Marker_mad8_adaptor : public Element_adaptor
{
public:
    Marker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Marker_mad8_adaptor();
};

class Drift_mad8_adaptor : public Element_adaptor
{
public:
    Drift_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Drift_mad8_adaptor();
};

class Sbend_mad8_adaptor : public Element_adaptor
{
public:
    Sbend_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Sbend_mad8_adaptor();
};

class Rbend_mad8_adaptor : public Element_adaptor
{
public:
    Rbend_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Rbend_mad8_adaptor();
};

class Quadrupole_mad8_adaptor : public Element_adaptor
{
public:
    Quadrupole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Quadrupole_mad8_adaptor();
};

class Sextupole_mad8_adaptor : public Element_adaptor
{
public:
    Sextupole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Sextupole_mad8_adaptor();
};

class Octupole_mad8_adaptor : public Element_adaptor
{
public:
    Octupole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Octupole_mad8_adaptor();
};

class Multipole_mad8_adaptor : public Element_adaptor
{
public:
    Multipole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Multipole_mad8_adaptor();
};

class Solenoid_mad8_adaptor : public Element_adaptor
{
public:
    Solenoid_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Solenoid_mad8_adaptor();
};

class Hkicker_mad8_adaptor : public Element_adaptor
{
public:
    Hkicker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Hkicker_mad8_adaptor();
};

class Vkicker_mad8_adaptor : public Element_adaptor
{
public:
    Vkicker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Vkicker_mad8_adaptor();
};

class Kicker_mad8_adaptor : public Element_adaptor
{
public:
    Kicker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Kicker_mad8_adaptor();
};

class Rfcavity_mad8_adaptor : public Element_adaptor
{
public:
    Rfcavity_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Rfcavity_mad8_adaptor();
};

class Elseparator_mad8_adaptor : public Element_adaptor
{
public:
    Elseparator_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Elseparator_mad8_adaptor();
};

class Hmonitor_mad8_adaptor : public Element_adaptor
{
public:
    Hmonitor_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Hmonitor_mad8_adaptor();
};

class Vmonitor_mad8_adaptor : public Element_adaptor
{
public:
    Vmonitor_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Vmonitor_mad8_adaptor();
};

class Monitor_mad8_adaptor : public Element_adaptor
{
public:
    Monitor_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Monitor_mad8_adaptor();
};

class Instrument_mad8_adaptor : public Element_adaptor
{
public:
    Instrument_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Instrument_mad8_adaptor();
};

class Ecollimator_mad8_adaptor : public Element_adaptor
{
public:
    Ecollimator_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Ecollimator_mad8_adaptor();
};

class Rcollimator_mad8_adaptor : public Element_adaptor
{
public:
    Rcollimator_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Rcollimator_mad8_adaptor();
};

#endif /* ELEMENT_ADAPTOR_H_ */
