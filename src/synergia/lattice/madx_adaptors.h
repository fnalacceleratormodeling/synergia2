#ifndef MADX_ADAPTORS_H_
#define MADX_ADAPTORS_H_

#include "synergia/lattice/element_adaptor.h"
#include "synergia/foundation/physical_constants.h"

class Marker_madx_adaptor : public Element_adaptor
{
public:
    Marker_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Marker_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Marker_madx_adaptor);

class Drift_madx_adaptor : public Element_adaptor
{
public:
    Drift_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Drift_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Drift_madx_adaptor);

class Sbend_madx_adaptor : public Element_adaptor
{
public:
    Sbend_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Sbend_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Sbend_madx_adaptor);

class Rbend_madx_adaptor : public Element_adaptor
{
public:
    Rbend_madx_adaptor();
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
    ~Rbend_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Rbend_madx_adaptor);

class Quadrupole_madx_adaptor : public Element_adaptor
{
public:
    static const char yoshida_propagator[];
    static const char basic_propagator[];
    Quadrupole_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Quadrupole_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Quadrupole_madx_adaptor);

class Sextupole_madx_adaptor : public Element_adaptor
{
public:
    Sextupole_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Sextupole_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Sextupole_madx_adaptor);

class Octupole_madx_adaptor : public Element_adaptor
{
public:
    Octupole_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Octupole_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Octupole_madx_adaptor);

class Multipole_madx_adaptor : public Element_adaptor
{
public:
    Multipole_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Multipole_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Multipole_madx_adaptor);

// thinpoles are an CHEF addon not found in MADX
class Thinpole_madx_adaptor : public Element_adaptor
{
public:
    Thinpole_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Thinpole_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Thinpole_madx_adaptor);

class Solenoid_madx_adaptor : public Element_adaptor
{
public:
    Solenoid_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Solenoid_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Solenoid_madx_adaptor);

class Hkicker_madx_adaptor : public Element_adaptor
{
public:
    Hkicker_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Hkicker_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Hkicker_madx_adaptor);

class Vkicker_madx_adaptor : public Element_adaptor
{
public:
    Vkicker_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Vkicker_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Vkicker_madx_adaptor);

class Kicker_madx_adaptor : public Element_adaptor
{
public:
    Kicker_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Kicker_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Kicker_madx_adaptor);

class Rfcavity_madx_adaptor : public Element_adaptor
{
public:
    Rfcavity_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    virtual void
    set_defaults(Lattice_element & lattice_element);
#if 0 // let CHEF set frequency
    virtual void
    set_derived_attributes_external(Lattice_element &lattice_element,
    		double lattice_length, double beta);
#endif // let CHEF set frequency
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Rfcavity_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Rfcavity_madx_adaptor);

class Elseparator_madx_adaptor : public Element_adaptor
{
public:
    Elseparator_madx_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Elseparator_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Elseparator_madx_adaptor);

class Hmonitor_madx_adaptor : public Element_adaptor
{
public:
    Hmonitor_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Hmonitor_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Hmonitor_madx_adaptor);

class Vmonitor_madx_adaptor : public Element_adaptor
{
public:
    Vmonitor_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Vmonitor_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Vmonitor_madx_adaptor);

class Monitor_madx_adaptor : public Element_adaptor
{
public:
    Monitor_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Monitor_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Monitor_madx_adaptor);

class Instrument_madx_adaptor : public Element_adaptor
{
public:
    Instrument_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Instrument_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Instrument_madx_adaptor);

class Ecollimator_madx_adaptor : public Element_adaptor
{
public:
    Ecollimator_madx_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Ecollimator_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Ecollimator_madx_adaptor);

class Rcollimator_madx_adaptor : public Element_adaptor
{
public:
    Rcollimator_madx_adaptor();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Rcollimator_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Rcollimator_madx_adaptor);

// Septum is an CHEF addon not found in MADX
class Septum_madx_adaptor : public Element_adaptor
{
public:
    Septum_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Septum_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Septum_madx_adaptor);

// Lambertson is an CHEF addon not found in MADX
class Lambertson_madx_adaptor : public Element_adaptor
{
public:
    Lambertson_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Lambertson_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Lambertson_madx_adaptor);


class Srot_madx_adaptor : public Element_adaptor
{
public:
    Srot_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Srot_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Srot_madx_adaptor);

class Dipedge_madx_adaptor : public Element_adaptor
{
public:
    Dipedge_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
    void
    serialize(Archive & ar, const unsigned int version);
    virtual
    ~Dipedge_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Dipedge_madx_adaptor);

class Nonlinearlens_madx_adaptor : public Element_adaptor
{
public:
    Nonlinearlens_madx_adaptor();
    Chef_elements
    get_chef_elements(Lattice_element const & lattice_element, double brho);
    template<class Archive>
    void
    serialize(Archive & ar, const unsigned int version);
    virtual
    ~Nonlinearlens_madx_adaptor();
};
BOOST_CLASS_EXPORT_KEY(Nonlinearlens_madx_adaptor);

#endif /* MADX_ADAPTORS_H_ */
