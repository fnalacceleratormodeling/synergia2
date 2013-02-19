#ifndef INDEPENDENT_STEPPER_ELEMENTS_H_
#define INDEPENDENT_STEPPER_ELEMENTS_H_
#include "synergia/simulation/stepper.h"

/// The Independent_stepper_elements class generates a constant number of
/// Independent_operator steps per thick element. Thin elements are assigned
/// a single step each. No collective effects are included.
class Independent_stepper_elements : public Stepper
{
private:
    void
    construct(int steps_per_element);
public:
    /// Construct an Independent_stepper
    /// @param lattice_sptr the Lattice
    /// @param map_order order for Chef_map operations
    /// @param steps_per_element the number of steps per thick element
    Independent_stepper_elements(Lattice_sptr lattice_sptr, int map_order,
            int steps_per_element);

    /// Deprecated. Construct an Independent_stepper
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param steps_per_element the number of steps per thick element
    Independent_stepper_elements(Lattice_simulator const& lattice_simulator,
            int steps_per_element);

    /// Default constructor for serialization use only
    Independent_stepper_elements();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Independent_stepper_elements();
};
BOOST_CLASS_EXPORT_KEY(Independent_stepper_elements);
typedef boost::shared_ptr<Independent_stepper_elements >
        Independent_stepper_elements_sptr;

#endif /* INDEPENDENT_STEPPER_ELEMENTS_H_ */
