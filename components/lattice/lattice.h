#ifndef LATTICE_H_
#define LATTICE_H_

#include <string>
#include <list>

#include "components/lattice/lattice_element.h"
#include "components/lattice/element_adaptor.h"
#include "components/foundation/reference_particle.h"
#include <boost/shared_ptr.hpp>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/version.hpp>

class Lattice
{
private:
    std::string name;
    Reference_particle *reference_particle_ptr;
    bool reference_particle_allocated;
    Lattice_elements elements;

public:
    Lattice();
    Lattice(std::string const& name);
    std::string const&
    get_name() const;
    void
    set_reference_particle(Reference_particle const& reference_particle);
    bool
    has_reference_particle() const;
    Reference_particle const&
    get_reference_particle() const;
    void
    append(Lattice_element const& element);
    void
    set_default_attributes(Element_adaptor_map const& element_adaptor_map);
    void
    set_default_attributes();
    Lattice_elements &
    get_elements();
    double
    get_length() const;
    double
    get_total_angle() const;
    void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(name)
                    & BOOST_SERIALIZATION_NVP(reference_particle_allocated);
            if (reference_particle_allocated) {
                ar & BOOST_SERIALIZATION_NVP(reference_particle_ptr);
            }
            ar & BOOST_SERIALIZATION_NVP(elements);
        }
    ~Lattice();
};

typedef boost::shared_ptr<Lattice > Lattice_sptr;

#endif /* LATTICE_H_ */
