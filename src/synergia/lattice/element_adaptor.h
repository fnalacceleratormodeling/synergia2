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

typedef boost::shared_ptr<Element_adaptor > Element_adaptor_sptr;

class Marker_mad8_adaptor : public Element_adaptor
{
public:
    Marker_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Marker_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Marker_mad8_adaptor);

class Drift_mad8_adaptor : public Element_adaptor
{
public:
    Drift_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Drift_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Drift_mad8_adaptor);

class Sbend_mad8_adaptor : public Element_adaptor
{
public:
    Sbend_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Sbend_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Sbend_mad8_adaptor);

class Rbend_mad8_adaptor : public Element_adaptor
{
public:
    Rbend_mad8_adaptor();
    virtual void
    set_defaults(Lattice_element & lattice_element);
    virtual void
    set_derived_attributes_internal(Lattice_element & lattice_element);
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Rbend_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Rbend_mad8_adaptor);

class Quadrupole_mad8_adaptor : public Element_adaptor
{
public:
    Quadrupole_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Quadrupole_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Quadrupole_mad8_adaptor);

class Sextupole_mad8_adaptor : public Element_adaptor
{
public:
    Sextupole_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Sextupole_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Sextupole_mad8_adaptor);

class Octupole_mad8_adaptor : public Element_adaptor
{
public:
    Octupole_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Octupole_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Octupole_mad8_adaptor);

class Multipole_mad8_adaptor : public Element_adaptor
{
public:
    Multipole_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Multipole_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Multipole_mad8_adaptor);

// thinpoles are an CHEF addon not found in MAD8
class Thinpole_mad8_adaptor : public Element_adaptor
{
public:
    Thinpole_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Thinpole_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Thinpole_mad8_adaptor);

class Solenoid_mad8_adaptor : public Element_adaptor
{
public:
    Solenoid_mad8_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Solenoid_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Solenoid_mad8_adaptor);

class Hkicker_mad8_adaptor : public Element_adaptor
{
public:
    Hkicker_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Hkicker_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Hkicker_mad8_adaptor);

class Vkicker_mad8_adaptor : public Element_adaptor
{
public:
    Vkicker_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Vkicker_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Vkicker_mad8_adaptor);

class Kicker_mad8_adaptor : public Element_adaptor
{
public:
    Kicker_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Kicker_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Kicker_mad8_adaptor);

class Rfcavity_mad8_adaptor : public Element_adaptor
{
public:
    Rfcavity_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Rfcavity_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Rfcavity_mad8_adaptor);

class Elseparator_mad8_adaptor : public Element_adaptor
{
public:
    Elseparator_mad8_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Elseparator_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Elseparator_mad8_adaptor);

class Hmonitor_mad8_adaptor : public Element_adaptor
{
public:
    Hmonitor_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Hmonitor_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Hmonitor_mad8_adaptor);

class Vmonitor_mad8_adaptor : public Element_adaptor
{
public:
    Vmonitor_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Vmonitor_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Vmonitor_mad8_adaptor);

class Monitor_mad8_adaptor : public Element_adaptor
{
public:
    Monitor_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Monitor_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Monitor_mad8_adaptor);

class Instrument_mad8_adaptor : public Element_adaptor
{
public:
    Instrument_mad8_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Instrument_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Instrument_mad8_adaptor);

class Ecollimator_mad8_adaptor : public Element_adaptor
{
public:
    Ecollimator_mad8_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Ecollimator_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Ecollimator_mad8_adaptor);

class Rcollimator_mad8_adaptor : public Element_adaptor
{
public:
    Rcollimator_mad8_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Rcollimator_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Rcollimator_mad8_adaptor);

// Septum is an CHEF addon not found in MAD8
class Septum_mad8_adaptor : public Element_adaptor
{
public:
    Septum_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Septum_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Septum_mad8_adaptor);

// Lambertson is an CHEF addon not found in MAD8
class Lambertson_mad8_adaptor : public Element_adaptor
{
public:
    Lambertson_mad8_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Lambertson_mad8_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Lambertson_mad8_adaptor);

#endif /* ELEMENT_ADAPTOR_H_ */
