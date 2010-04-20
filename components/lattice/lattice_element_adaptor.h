#ifndef LATTICE_ELEMENT_ADAPTOR_H_
#define LATTICE_ELEMENT_ADAPTOR_H_

#include "components/lattice/lattice_element.h"
#include <boost/shared_ptr.hpp>

class Lattice_element_adaptor
{
public:
    Lattice_element_adaptor();
    virtual void
    set_default_attributes(Lattice_element & lattice_element) = 0;
    virtual
    ~Lattice_element_adaptor();
};

typedef boost::shared_ptr<Lattice_element_adaptor >
        Lattice_element_adaptor_sptr;

class Marker_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Marker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Marker_mad8_adaptor();
};

class Drift_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Drift_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Drift_mad8_adaptor();
};

class Sbend_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Sbend_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Sbend_mad8_adaptor();
};

class Rbend_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Rbend_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Rbend_mad8_adaptor();
};

class Quadrupole_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Quadrupole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Quadrupole_mad8_adaptor();
};

class Sextupole_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Sextupole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Sextupole_mad8_adaptor();
};

class Octupole_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Octupole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Octupole_mad8_adaptor();
};

class Multipole_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Multipole_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Multipole_mad8_adaptor();
};

class Solenoid_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Solenoid_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Solenoid_mad8_adaptor();
};

class Hkicker_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Hkicker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Hkicker_mad8_adaptor();
};

class Vkicker_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Vkicker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Vkicker_mad8_adaptor();
};

class Kicker_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Kicker_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Kicker_mad8_adaptor();
};

class Rfcavity_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Rfcavity_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Rfcavity_mad8_adaptor();
};

class Elseparator_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Elseparator_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Elseparator_mad8_adaptor();
};

class Hmonitor_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Hmonitor_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Hmonitor_mad8_adaptor();
};

class Vmonitor_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Vmonitor_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Vmonitor_mad8_adaptor();
};

class Monitor_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Monitor_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Monitor_mad8_adaptor();
};

class Instrument_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Instrument_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Instrument_mad8_adaptor();
};

class Ecollimator_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Ecollimator_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Ecollimator_mad8_adaptor();
};

class Rcollimator_mad8_adaptor : public Lattice_element_adaptor
{
public:
    Rcollimator_mad8_adaptor();
    virtual void
    set_default_atributes(Lattice_element & lattice_element);
    virtual
    ~Rcollimator_mad8_adaptor();
};

#endif /* LATTICE_ELEMENT_ADAPTOR_H_ */
