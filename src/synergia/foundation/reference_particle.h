#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_serialization.h"
#include "synergia/foundation/four_momentum.h"

/// Reference_particle stores the four momentum of the reference frame
/// with respect to  the lab frame (defined to be along the axis of the
/// accelerator) as well as the six-dimensional state vector of the the
/// reference particle in the reference frame. Reference particle
/// also keeps track of the total path length of the reference particle
/// trajectory.
class Reference_particle
{
private:
    int charge;
    Four_momentum four_momentum;
    MArray1d state;
    int repetition;
    double repetition_length;
    double s_n;
public:
    /// Default constructor for internal use only
    Reference_particle();

    /// Construct a Reference_particle with a given mass and total energy.
    /// @param mass in GeV/c^2
    /// @param charge in units of e
    /// @param total_energy in GeV in the lab frame
    Reference_particle(int charge, double mass, double total_energy);

    /// Construct a Reference_particle with a given four momentum.
    /// @param charge in units of e
    /// @param four_momentum in the lab frame
    Reference_particle(int charge, Four_momentum const& four_momentum);

    /// Construct a Reference_particle with a given four momentum and state
    /// in the reference frame.
    /// @param charge in units of e
    /// @param four_momentum in the lab frame
    /// @param state is a six-dimensional state vector
    Reference_particle(int charge, Four_momentum const& four_momentum,
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

    /// Increment the trajectory length.
    /// @param length in m
    void
    increment_trajectory(double length);

    /// Start a new repetition
    void
    start_repetition();

    /// Manually set trajectory parameters
    /// @param repetition starting at 0
    /// @param repetition_length in m
    /// @param s in m
    void
    set_trajectory(int repetition, double repetition_length, double s);

    /// Return the Reference_particle charge in units of e
    int
    get_charge() const;

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

    /// Get the total path length in m of the reference
    /// particle trajectory
    double
    get_trajectory_length() const;

    /// Get the distance traveled in m since the beginning
    /// of the current repetition.
    double
    get_s_n() const;

    /// Get the number of repetition.
    int
    get_repetition() const;

    /// Get the repetition length in m.
    double
    get_repetition_length() const;

    /// Check equality to the given tolerance
    /// @param reference_particle another Reference_particle
    /// @param tolerance fractional accuracy
    bool
    equal(Reference_particle const& reference_particle, double tolerance) const;

    /// Serialization support
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

#endif /* REFERENCE_PARTICLE_H_ */
