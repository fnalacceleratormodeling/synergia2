#ifndef LATTICE_H_
#define LATTICE_H_

#include <string>
#include <list>

#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/element_adaptor.h"
#include "synergia/foundation/reference_particle.h"
#include <boost/shared_ptr.hpp>

#include "synergia/utils/serialization.h"

/// The Lattice class contains an abstract representation of an ordered
/// set of Lattice_elements. Each element of the Lattice is unique.
class Lattice
{
private:
    std::string name;
    Reference_particle *reference_particle_ptr;
    bool reference_particle_allocated;
    Lattice_elements elements;

public:
    /// Construct a Lattice object without a name.
    Lattice();

    /// Construct a Lattice object with a name
    /// @param name an arbitrary name
    Lattice(std::string const& name);

    /// Get the Lattice name
    std::string const&
    get_name() const;

    /// Set the Lattice reference particle
    /// @param reference_particle a Reference_particle
    void
    set_reference_particle(Reference_particle const& reference_particle);

    /// Determine whether the Lattice has a reference particle
    bool
    has_reference_particle() const;

    /// Get the Lattice reference particle
    Reference_particle const&
    get_reference_particle() const;

    /// Append a copy of a Lattice_element.
    /// @param element a Lattice_element
    void
    append(Lattice_element const& element);

    /// Set the default Element_adaptor_map to be used for elements in the
    /// Lattice
    /// @param element_adaptor_map an Element_adaptor_map
    void
    set_default_attributes(Element_adaptor_map const& element_adaptor_map);

    /// Derive internal attributes where necessary
    void
    derive_internal_attributes(Element_adaptor_map const& element_adaptor_map);

    /// Derive external attributes where necessary
    void
    derive_external_attributes(Element_adaptor_map const& element_adaptor_map);

    /// Complete all attribute updates. Includes defaults and derivations.
    void
    complete_attributes(Element_adaptor_map const& element_adaptor_map);

    /// Get the list of elements in the Lattice
    Lattice_elements &
    get_elements();

    /// Get the combined length of all the elements in the Lattice
    double
    get_length() const;

    /// Get the total angle in radians subtended by all the elements in the
    /// Lattice
    double
    get_total_angle() const;

    /// Print a human-readable summary of the elements in the Lattice.
    /// The Python version of this function is named "print_".
    void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~Lattice();
};

typedef boost::shared_ptr<Lattice > Lattice_sptr;

#endif /* LATTICE_H_ */
