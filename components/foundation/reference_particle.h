#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

#include "utils/multi_array_typedefs.h"
#include "components/foundation/four_momentum.h"

/// Reference_particle stores the four momentum of the reference frame
/// with respect to  the lab frame (defined to be along the axis of the
/// accelerator) as well as the six-dimensional state vector of the the
/// reference particle in the reference frame.
class Reference_particle
{
private:
    Four_momentum four_momentum;
    MArray1d state;
public:
    /// Construct a Reference_particle with a given mass and total energy.
    /// @param mass in GeV/c^2
    /// @param total_energy in GeV in the lab frame
    Reference_particle(double mass, double total_energy);

    /// Construct a Reference_particle with a given four momentum.
    /// @param four_momentum in the lab frame
    Reference_particle(Four_momentum const& four_momentum);

    /// Construct a Reference_particle with a given four momentum and state
    /// in the reference frame.
    /// @param four_momentum in the lab frame
    /// @param state is a six-dimensional state vector
    Reference_particle(Four_momentum const& four_momentum,
            Const_MArray1d_ref state);

    /// Set the four momentum.
    /// @param four_momentum in the lab frame
    void
    set_four_momentum(Four_momentum const& four_momentum);

    /// Set the state vector in the reference frame.
    /// @param state is a six-dimensional state vector
    void
    set_state(Const_MArray1d_ref state);

    /// Set the total energy.
    /// @param total_energy in GeV in the lab frame
    void
    set_total_energy(double total_energy);

    /// Get the four momentum in the lab frame.
    Four_momentum const &
    get_four_momentum() const;

    /// Get the six-dimensional state vector in the reference frame.
    Const_MArray1d_ref
    get_state() const;

    /// Get the relativistic beta in the lab frame.
    double
    get_beta() const;

    /// Get the relativistic gamma in the lab frame.
    double
    get_gamma() const;

    /// Get the momentum in GeV/c in the lab frame.
    double
    get_momentum() const;

    /// Get the total energy in GeV in the lab frame.
    double
    get_total_energy() const;

    /// Check equality to the given tolerance
    /// @param reference_particle another Reference_particle
    /// @param tolerance fractional accuracy
    bool
    equal(Reference_particle const& reference_particle, double tolerance) const;
};

#endif /* REFERENCE_PARTICLE_H_ */
