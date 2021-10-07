#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element_processor.h"

//#include "mad8_adaptor_map.h"
//#include "madx_adaptor_map.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

Lattice::Lattice() 
    : name("")
    , reference_particle()
    , elements()
    , updated { true, true, true }
{
}

Lattice::Lattice(std::string const & name) 
    : name(name)
    , reference_particle()
    , elements()
    , updated { true, true, true }
{
}

Lattice::Lattice(Lattice const & o) 
    : name(o.name)
    , reference_particle(o.reference_particle)
    , elements(o.elements)
    , updated(o.updated)
{
    for (auto & e : elements) e.set_lattice(*this);
}

Lattice::Lattice(Lattice && o) noexcept
    : name(std::move(o.name))
    , reference_particle(std::move(o.reference_particle))
    , elements(std::move(o.elements))
    , updated(std::move(o.updated))
{
    for (auto & e : elements) e.set_lattice(*this);
}

Lattice&
Lattice::operator=(Lattice const& o)
{
    name = o.name;
    reference_particle = o.reference_particle;
    elements = o.elements;
    updated = o.updated;

    for (auto & e : elements) e.set_lattice(*this);

    return *this;
}


Lattice::Lattice(Lsexpr const & lsexpr) 
    : name("")
    , reference_particle()
    , elements()
    , updated { true, true, true }
{
    for (auto const& lse : lsexpr)
    {
        if (lse.is_labeled()) 
        {
            if (lse.get_label() == "name") 
            {
                name = lse.get_string();
            } 
#if 0
            else if (it->get_label() == "type") 
            {
                std::string lctype(it->get_string());
                std::transform(lctype.begin(), lctype.end(), lctype.begin(), ::tolower);
                if (lctype == "mad8") 
                {
                    element_adaptor_map_sptr = boost::shared_ptr<Element_adaptor_map>(new Mad8_adaptor_map);
                } 
                else if (lctype == "madx") 
                {
                    element_adaptor_map_sptr = boost::shared_ptr<Element_adaptor_map>(new MadX_adaptor_map);
                } 
                else 
                {
                    throw std::runtime_error("Lattice: adaptor map type " +
                                             it->get_string() + " not handled");
                }
            } 
#endif
            else if (lse.get_label() == "reference_particle") 
            {
                reference_particle = Reference_particle(lse);
            } 
            else if (lse.get_label() == "elements") 
            {
                for (auto const& ele : lse)
                {
                    append(Lattice_element(ele));
                }
            }
        }
    }

#if 0
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
                reference_particle = new Reference_particle(*it);
                reference_particle_allocated = true;
            } else if (it->get_label() == "elements") {
                for (Lsexpr::const_iterator_t eit = it->begin();
                     eit != it->end(); ++eit) {
                    append(Lattice_element(*eit));
                }
            }
        }
    }
#endif
}

#if 0
Lsexpr
Lattice::as_lsexpr() const
{
    Lsexpr retval;
#if 0
    retval.push_back(Lsexpr(name, "name"));
    retval.push_back(Lsexpr(element_adaptor_map_sptr->get_label(),
                            "type"));
    if (reference_particle) {
        Lsexpr ref_lsexpr(reference_particle->as_lsexpr());
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
#endif
    return retval;
}
#endif

Lattice::update_flags_t
Lattice::update()
{
    update_flags_t res = updated;

    if ( updated.ref || updated.structure || updated.element)
    {
        // TODO: reconstruct the chef lattice
        // ...

        updated = { false, false, false };
    }

    return res;
}

void
Lattice::append(Lattice_element const & element)
{
    elements.push_back(Lattice_element_processor::process(element));
    elements.back().set_lattice(*this);
    updated.structure = true;
}

#if 0
void
Lattice::derive_external_attributes()
{
#if 0
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
        double beta = reference_particle->get_beta();
        double lattice_length = get_length();
        for (Lattice_elements::const_iterator it = elements.begin();
                it != elements.end(); ++it) {
            if ((*it)->get_needs_external_derive()) {
                element_adaptor_map_sptr->get_adaptor((*it)->get_type())->set_derived_attributes_external(
                        *(*it), lattice_length, beta);
            }
        }
    }
#endif
}
#endif

void
Lattice::set_all_double_attribute(
        std::string const & name, 
        double value,
        bool increment_revision)
{
    for (auto & e : elements)
        e.set_double_attribute(name, value, increment_revision);

    updated.element = true;
}

void
Lattice::set_all_string_attribute(
        std::string const & name,
        std::string const & value, 
        bool increment_revision )
{
    for (auto & e : elements)
        e.set_string_attribute(name, value, increment_revision);

    updated.element = true;
}

std::list<Lattice_element> const &
Lattice::get_elements() const
{
    return elements;
}

std::list<Lattice_element> &
Lattice::get_elements()
{
    return elements;
}

void
Lattice::reset_all_markers()
{
    for(auto& e : elements) 
        e.reset_markers();
}

double
Lattice::get_length() const
{
    double length = 0.0;
    for (auto const & e : elements) length += e.get_length();
    return length;
}

double
Lattice::get_total_angle() const
{
    double angle = 0.0;
    for (auto const & e : elements) angle += e.get_bend_angle();
    return angle;
}

std::string
Lattice::as_string() const
{
    std::stringstream sstream;
    sstream << name << ":\n";
    for (auto const & e : elements) sstream << e.as_string() << std::endl;
    return sstream.str();
}

void
Lattice::print(Logger & logger) const
{
    logger(LoggerV::DEBUG) << as_string() << std::endl;
}

