#ifndef INDEPENDENT_STEPPER_H_
#define INDEPENDENT_STEPPER_H_
#include "synergia/simulation/stepper.h"

/// The Independent_stepper class generates evenly-spaced Independent_operator
/// steps through a Lattice. No collective effects are included.
class Independent_stepper : public Stepper
{
private:
    void
    construct(int num_steps);
public:
    /// Construct an Independent_stepper
    /// @param lattice_sptr the Lattice
    /// @param map_order order for Chef_map operations
    /// @param num_steps the number of steps to take in the Lattice
    Independent_stepper(Lattice_sptr lattice_sptr, int map_order,
            int num_steps);

    /// Deprecated. Construct an Independent_stepper
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param num_steps the number of steps to take in the Lattice
    Independent_stepper(Lattice_simulator const& lattice_simulator,
            int num_steps);

    /// Default constructor for serialization use only
    Independent_stepper();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Independent_stepper();

};
BOOST_CLASS_EXPORT_KEY(Independent_stepper);
typedef boost::shared_ptr<Independent_stepper > Independent_stepper_sptr;

#endif /* INDEPENDENT_STEPPER_H_ */
