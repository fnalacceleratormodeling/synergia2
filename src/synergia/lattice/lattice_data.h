#ifndef LATTICE_DATA_H_
#define LATTICE_DATA_H_

#include <string>
#include <list>

#include "synergia/lattice/lattice_element.h"
#include "synergia/foundation/reference_particle.h"

//#include "synergia/lattice/diagnostics_apertures_loss.h"
//#include "synergia/lattice/element_adaptor.h"
//#include "synergia/lattice/element_adaptor_map.h"

#include "synergia/utils/cereal.h"

/// The Lattice_data class contains an abstract representation of an ordered
/// set of objects of type Lattice_element.
/// Each element of the Lattice_data is unique.
class Lattice_data
{

private:

    std::string name;

    Reference_particle reference_particle;
    Lattice_elements elements;

    bool dirty;

    //Element_adaptor_map_sptr element_adaptor_map_sptr;

    //Diagnostics_losses diagnostics_loss_list;
    //bool have_loss_diagnostics;

public:

    /// Construct a Lattice_data object without a name.
    /// Defaults to interpreting elements as Mad8 elements
    Lattice_data();

    /// Copies of Lattice_datas contain copies of elements
    Lattice_data(Lattice_data const & lattice);

    /// Construct a Lattice_data object with a name
    /// Defaults to interpreting elements as Mad8 elements
    /// @param name an arbitrary name
    explicit Lattice_data(std::string const & name);

    /// Construct a Lattice_data from the Lsexpr representation
    /// @param lsexpr representation
    explicit Lattice_data(Lsexpr const & lsexpr);

    /// Extract an Lsexpr representation of the Lattice_data
    Lsexpr
    as_lsexpr() const;

    /// Get the Lattice_data name
    std::string const &
    get_name() const;

    /// Set the Lattice_data reference particle
    /// @param reference_particle a Reference_particle
    void
    set_reference_particle(Reference_particle const & reference_particle);

    /// Get the Lattice_data reference particle (const)
    Reference_particle const &
    get_reference_particle() const;

    /// Append a copy of a Lattice_element.
    /// @param element a Lattice_element
    void
    append(Lattice_element const & element);

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
    set_all_double_attribute(
            std::string const& name, double value,
            bool increment_revision = true);

    /// Set the value of the named string attribute on all elements
    /// @param name attribute name
    /// @param value attribute value
    /// @param increment_revision can be set to false for attributes that do not affect dynamics
    void
    set_all_string_attribute(
            std::string const & name, std::string const & value,
            bool increment_revision = true);

    /// Get the list of elements in the Lattice_data
    Lattice_elements const &
    get_elements() const;

    /// Get the combined length of all the elements in the Lattice_data
    double
    get_length() const;

    /// Get the total angle in radians subtended by all the elements in the
    /// Lattice_data
    double
    get_total_angle() const;

    /// Return a human-readable summary of the elements in the Lattice_data.
    std::string
    as_string() const;

    /// Print a human-readable summary of the elements in the Lattice_data.
    /// The Python version of this function is named "print_".
    void
    print() const;

    template<class Archive>
    void
    serialize(Archive & ar, const unsigned int version);
};

typedef std::shared_ptr<Lattice_data> Lattice_data_sptr; // syndoc:include

#endif /* LATTICE_H_ */
