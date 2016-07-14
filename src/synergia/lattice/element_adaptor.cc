#include <stdexcept>
#include <sstream>
#include "element_adaptor.h"
#include "synergia/foundation/math_constants.h"

#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <beamline/beamline_elements.h>
#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif

Element_adaptor::Element_adaptor() :
        default_element_sptr(new Lattice_element)
{
}

Lattice_element_sptr
Element_adaptor::get_default_element_sptr()
{
    return default_element_sptr;
}

Lattice_element &
Element_adaptor::get_default_element()
{
    return *default_element_sptr;
}

void
Element_adaptor::set_double_default(Lattice_element & lattice_element,
        std::string const& name, double value)
{
    if (!lattice_element.has_double_attribute(name)) {
        lattice_element.set_double_attribute(name, value);
    }
}

void
Element_adaptor::set_string_default(Lattice_element & lattice_element,
        std::string const& name, std::string const& value)
{
    if (!lattice_element.has_string_attribute(name)) {
        lattice_element.set_string_attribute(name, value);
    }
}

void
Element_adaptor::set_defaults(Lattice_element & lattice_element)
{
    lattice_element.set_default_element(get_default_element_sptr());
}

void
Element_adaptor::set_derived_attributes_internal(
        Lattice_element & lattice_element)
{
}

void
Element_adaptor::set_derived_attributes_external(
        Lattice_element & lattice_element, double lattice_length, double beta)
{
}

Chef_elements
Element_adaptor::get_chef_elements(Lattice_element const& lattice_element,
        double brho)
{
    throw(runtime_error(
            "Element_adaptor: " + lattice_element.get_type() + " not handled"));
}

template<class Archive>
    void
    Element_adaptor::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(default_element_sptr);
    }

template
void
Element_adaptor::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Element_adaptor::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Element_adaptor::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Element_adaptor::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Element_adaptor::~Element_adaptor()
{
}
