#ifndef LATTICE_H_
#define LATTICE_H_

#include <string>
#include <list>

#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/element_adaptor.h"
#include "synergia/lattice/element_adaptor_map.h"
#include "synergia/lattice/diagnostics_apertures_loss.h"
#include "synergia/foundation/reference_particle.h"
#include <boost/shared_ptr.hpp>

#include "synergia/utils/serialization.h"

/// The Lattice class contains an abstract representation of an ordered
/// set of objects of type Lattice_element.
/// Each element of the Lattice is unique.
class Lattice
{
private:
    std::string name;
    bool reference_particle_allocated;
    Reference_particle *reference_particle_ptr;
    Lattice_elements elements;
    Element_adaptor_map_sptr element_adaptor_map_sptr;
    Diagnostics_losses diagnostics_loss_list;
    bool have_loss_diagnostics;
   

public:
    /// Construct a Lattice object without a name.
    /// Defaults to interpreting elements as Mad8 elements
    Lattice();

    /// Construct a Lattice object with a name
    /// Defaults to interpreting elements as Mad8 elements
    /// @param name an arbitrary name
    Lattice(std::string const& name);

    /// Construct a Lattice object with a name
    /// @param name an arbitrary name
    /// @param element_adaptor_map_sptr an Element_adaptor_map for interpreting elements
    Lattice(std::string const& name,
            Element_adaptor_map_sptr element_adaptor_map_sptr);

    /// Construct a Lattice from the Lsexpr representation
    /// @param lsexpr representation
    Lattice(Lsexpr const& lsexpr);

    /// Extract an Lsexpr representation of the Lattice
    Lsexpr
    as_lsexpr() const;

    /// Copies of Lattices contain copies of elements
     Lattice(Lattice const& lattice);

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
    Reference_particle &
    get_reference_particle();

    /// Get the Lattice reference particle (const)
    Reference_particle const &
    get_reference_particle() const;

    /// Append a copy of a Lattice_element.
    /// @param element a Lattice_element
    void
    append(Lattice_element const& element);

    /// Set the defaults in  elements of the Lattice
    void
    set_defaults();

    /// Derive internal attributes where necessary
    void
    derive_internal_attributes();

    /// Derive external attributes where necessary
    void
    derive_external_attributes();

    /// Complete all attribute updates. Includes defaults and derivations.
    void
    complete_attributes();

    /// Set the value of the named double attribute on all elements
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_all_double_attribute(std::string const& name, double value,
            bool increment_revision = true);

    /// Set the value of the named string attribute on all elements
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_all_string_attribute(std::string const& name, std::string const& value,
            bool increment_revision = true);

    /// Get the list of elements in the Lattice
    Lattice_elements &
    get_elements();

    /// Get the Element_adaptor_map
    Element_adaptor_map &
    get_element_adaptor_map();

    /// Get the Element_adaptor_map as a shared pointer
    Element_adaptor_map_sptr
    get_element_adaptor_map_sptr();
    Element_adaptor_map_sptr
    get_element_adaptor_map_sptr() const;

    /// Get the combined length of all the elements in the Lattice
    double
    get_length() const;

    /// Get the total angle in radians subtended by all the elements in the
    /// Lattice
    double
    get_total_angle() const;

    /// Return a human-readable summary of the elements in the Lattice.
    std::string
    as_string() const;

    bool 
    get_have_loss_diagnostics() const;
    
    Diagnostics_losses
    get_loss_diagnostics_list();
  
  //  Diagnostics_apertures_losses
  //  get_diagnostics_list() const;
    
    void
    add_loss_diagnostics(Diagnostics_loss_sptr diagnostics_sptr);
    
    /// Print a human-readable summary of the elements in the Lattice.
    /// The Python version of this function is named "print_".
    void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~Lattice();
};

typedef boost::shared_ptr<Lattice > Lattice_sptr; // syndoc:include

#endif /* LATTICE_H_ */
