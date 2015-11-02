#include "lattice.h"
#include "mad8_adaptor_map.h"
#include "madx_adaptor_map.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

Lattice::Lattice() :
        name(""), reference_particle_allocated(false), reference_particle_ptr(
                0), elements(), element_adaptor_map_sptr(new Mad8_adaptor_map)
{
}

Lattice::Lattice(std::string const& name) :
        name(name), reference_particle_allocated(false), elements(), element_adaptor_map_sptr(
                new Mad8_adaptor_map)
{
}

Lattice::Lattice(std::string const& name,
        Element_adaptor_map_sptr element_adaptor_map_sptr) :
        name(name), reference_particle_allocated(false), elements(), element_adaptor_map_sptr(
                element_adaptor_map_sptr)
{
}

Lattice::Lattice(Lattice const& lattice) :
        name(lattice.name), reference_particle_allocated(false), elements()
{
    for (Lattice_elements::const_iterator it = lattice.elements.begin();
            it != lattice.elements.end(); ++it) {
        Lattice_element_sptr lattice_element_sptr(new Lattice_element(**it));
        elements.push_back(lattice_element_sptr);
        lattice_element_sptr->set_lattice(*this);
    }
    if (lattice.reference_particle_allocated) {
        reference_particle_ptr = new Reference_particle(
                *lattice.reference_particle_ptr);
        reference_particle_allocated = true;
    }
}

Lattice::Lattice(Lsexpr const& lsexpr) :
      name("")
    , reference_particle_allocated(false)
    , reference_particle_ptr(0)
    , elements()
    , element_adaptor_map_sptr()
{
    for (Lsexpr::const_iterator_t it = lsexpr.begin(); it != lsexpr.end();
         ++it) {
        if (it->is_labeled()) {
            if (it->get_label() == "name") {
                name = it->get_string();
            } else if (it->get_label() == "type") {
                std::string lctype(it->get_string());
                std::transform(lctype.begin(), lctype.end(), lctype.begin(),
                               ::tolower);
                if (lctype == "mad8") {
                    element_adaptor_map_sptr = boost::shared_ptr<Element_adaptor_map>(new Mad8_adaptor_map);
                } else if (lctype == "madx") {
                    element_adaptor_map_sptr = boost::shared_ptr<Element_adaptor_map>(new MadX_adaptor_map);
                } else {
                    throw std::runtime_error("Lattice: adaptor map type " +
                                             it->get_string() + " not handled");
                }
            } else if (it->get_label() == "reference_particle") {
                reference_particle_ptr = new Reference_particle(*it);
                reference_particle_allocated = true;
            } else if (it->get_label() == "elements") {
                for (Lsexpr::const_iterator_t eit = it->begin();
                     eit != it->end(); ++eit) {
                    append(Lattice_element(*eit));
                }
            }
        }
    }
}

Lsexpr
Lattice::as_lsexpr() const
{
    Lsexpr retval;
    Lsexpr name_lsexpr(name);
    name_lsexpr.set_label("name");
    retval.push_back(name_lsexpr);
    Lsexpr type_lsexpr(element_adaptor_map_sptr->get_label());
    type_lsexpr.set_label("type");
    retval.push_back(type_lsexpr);
    if (reference_particle_ptr) {
        Lsexpr ref_lsexpr(reference_particle_ptr->as_lsexpr());
        ref_lsexpr.set_label("reference_particle");
        retval.push_back(ref_lsexpr);
    }
    Lsexpr elements_lsexpr;
    elements_lsexpr.set_label("elements");
    for(Lattice_elements::const_iterator it = elements.begin();
        it != elements.end(); ++it) {
        Lsexpr element_lsexpr((*it)->as_lsexpr());
        elements_lsexpr.push_back(element_lsexpr);
    }
    retval.push_back(elements_lsexpr);
    return retval;
}

std::string const&
Lattice::get_name() const
{
    return name;
}

void
Lattice::set_reference_particle(Reference_particle const& reference_particle)
{
    reference_particle_ptr = new Reference_particle(reference_particle);
    reference_particle_allocated = true;
}

bool
Lattice::has_reference_particle() const
{
    return reference_particle_allocated;
}

Reference_particle &
Lattice::get_reference_particle()
{
    if (!reference_particle_allocated) {
        throw std::runtime_error("Lattice: no reference_particle available");
    }
    return *reference_particle_ptr;
}

void
Lattice::append(Lattice_element const& element)
{
    Lattice_element_sptr element_sptr(new Lattice_element(element));
    elements.push_back(element_sptr);
    element_sptr->set_lattice(*this);
}

void
Lattice::set_defaults()
{
    for (Lattice_elements::const_iterator it = elements.begin();
            it != elements.end(); ++it) {
        if (element_adaptor_map_sptr->has_adaptor((*it)->get_type())) {
            element_adaptor_map_sptr->get_adaptor((*it)->get_type())->set_defaults(
                    *(*it));
        }
    }
}

void
Lattice::derive_internal_attributes()
{
    for (Lattice_elements::const_iterator it = elements.begin();
            it != elements.end(); ++it) {
        if ((*it)->get_needs_internal_derive()) {
            element_adaptor_map_sptr->get_adaptor((*it)->get_type())->set_derived_attributes_internal(
                    *(*it));
        }
    }
}

void
Lattice::derive_external_attributes()
{
    bool needed = false;
    for (Lattice_elements::const_iterator it = elements.begin();
            it != elements.end(); ++it) {
        if ((*it)->get_needs_external_derive()) {
            needed = true;
        }
    }
    if (needed) {
        if (!reference_particle_allocated) {
            throw std::runtime_error(
                    "Lattice::derive_external_attributes requires a reference_particle");
        }
        double beta = reference_particle_ptr->get_beta();
        double lattice_length = get_length();
        for (Lattice_elements::const_iterator it = elements.begin();
                it != elements.end(); ++it) {
            if ((*it)->get_needs_external_derive()) {
                element_adaptor_map_sptr->get_adaptor((*it)->get_type())->set_derived_attributes_external(
                        *(*it), lattice_length, beta);
            }
        }
    }
}

void
Lattice::complete_attributes()
{
    set_defaults();
    derive_internal_attributes();
    derive_external_attributes();
}

void
Lattice::set_all_double_attribute(std::string const& name, double value,
        bool increment_revision)
{
    for (Lattice_elements::iterator it = elements.begin(); it != elements.end();
            ++it) {
        (*it)->set_double_attribute(name, value, increment_revision);
    }
}

void
Lattice::set_all_string_attribute(std::string const& name,
        std::string const& value, bool increment_revision)
{
    for (Lattice_elements::iterator it = elements.begin(); it != elements.end();
            ++it) {
        (*it)->set_string_attribute(name, value, increment_revision);
    }
}

Lattice_elements &
Lattice::get_elements()
{
    return elements;
}

Element_adaptor_map &
Lattice::get_element_adaptor_map()
{
    return *element_adaptor_map_sptr;
}

Element_adaptor_map_sptr
Lattice::get_element_adaptor_map_sptr()
{
    return element_adaptor_map_sptr;
}

double
Lattice::get_length() const
{
    double length = 0.0;
    for (Lattice_elements::const_iterator it = elements.begin();
            it != elements.end(); ++it) {
        length += (*it)->get_length();
    }
    return length;
}

double
Lattice::get_total_angle() const
{
    double angle = 0.0;
    for (Lattice_elements::const_iterator it = elements.begin();
            it != elements.end(); ++it) {
        angle += (*it)->get_bend_angle();
    }
    return angle;
}

std::string
Lattice::as_string() const
{
    std::stringstream sstream;
    sstream << name << ":\n";
    for (Lattice_elements::const_iterator it = elements.begin();
            it != elements.end(); ++it) {
        sstream << (*it)->as_string();
        sstream << std::endl;
    }
    return sstream.str();
}

void
Lattice::print() const
{
    std::cout << as_string() << std::endl;
}

Lattice::~Lattice()
{
    if (reference_particle_allocated) {
        delete reference_particle_ptr;
    }

}

template<class Archive>
    void
    Lattice::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(name)& BOOST_SERIALIZATION_NVP(reference_particle_allocated);
        if (reference_particle_allocated) {
            ar & BOOST_SERIALIZATION_NVP(reference_particle_ptr);
        }
        ar & BOOST_SERIALIZATION_NVP(elements);
        ar & BOOST_SERIALIZATION_NVP(element_adaptor_map_sptr);
    }

template
void
Lattice::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lattice::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Lattice::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lattice::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
