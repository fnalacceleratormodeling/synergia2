#ifndef BUNCH_H_
#define BUNCH_H_

#include <mpi.h>
#include "utils/multi_array_typedefs.h"
#include "components/foundation/reference_particle.h"
#include "utils/commxx.h"
#include "components/bunch/fixed_t_z_converter.h"

/// Bunch represents a macroparticle bunch distributed across the processors
/// in a communicator.
class Bunch
{
public:
    enum State
    {
        fixed_z = 1, fixed_t = 2
    };
    static const int x = 0;
    static const int xp = 1;
    static const int y = 2;
    static const int yp = 3;
    static const int z = 4;
    static const int zp = 5;
    static const int t = 4;
    static const int tp = 5;
    static const int id = 6;
private:
    Reference_particle reference_particle;
    int particle_charge;
    MArray2d *local_particles;
    int local_num, total_num;
    double real_num;
    State state;
    Commxx comm;
    Fixed_t_z_converter *converter_ptr;
    Fixed_t_z_zeroth default_converter;
    void
    assign_ids(int local_offset);
    void
    construct(int particle_charge, int total_num, double real_num);
public:
    /// Construct a bunch. Allocates memory for the particles and assigns
    /// partice ID's, but does not fill the phase space values in any way.
    /// @param reference_particle the reference particle for the bunch.
    /// @param total_num the total number of macroparticles in the bunch
    /// @param real_num the number of real particles represented by the bunch.
    /// @param comm the communicator.
    Bunch(Reference_particle const& reference_particle, int total_num,
            double real_num, Commxx const& comm);

    /// Construct a bunch. Allocates memory for the particles and assigns
    /// partice ID's, but does not fill the phase space values in any way.
    /// @param reference_particle the reference particle for the bunch.
    /// @param total_num the total number of macroparticles in the bunch
    /// @param real_num the number of real particles represented by the bunch.
    /// @param comm the communicator.
    /// @param particle_charge in units of e.
    Bunch(Reference_particle const& reference_particle, int total_num,
            double real_num, Commxx const& comm, int particle_charge);

    Bunch(Bunch const& bunch);
    Bunch &
    operator=(Bunch const& bunch);

    /// Set the particle charge
    /// @param particle_charge in units of e.
    void
    set_particle_charge(int particle_charge);

    /// Set the number of real particles represented by the bunch.
    /// @param real_num the new real number of particles
    void
    set_real_num(double real_num);

    /// Reduce (set) the number of particles on this processor. The number
    /// of particles can only be lowered by this member function. (In order
    /// to add new particles, create another Bunch and use the inject member.)
    /// The total number and real number for the bunch will not be correct
    /// until update_total_num() is called. The real number will scale to
    /// reflect the change in the total number. n.b.: The only way to change
    /// the total number after the bunch has been created is to change the
    /// local numbers on each processor.
    /// @param local_num the new number of particles on this processor
    void
    set_local_num(int local_num);

    /// Update the total number and real number of particles after the local
    /// number has been changed. Requires communication.
    void
    update_total_num();

    /// Set the Fixed_t_z_converter class to be used when converting between
    /// fixed-t and fixed-z representations.
    /// @param converter the converter class
    void
    set_converter(Fixed_t_z_converter &converter);

    /// Convert to (fixed-t or fixed-z) state if necessary. Does nothing if
    /// the bunch is already in the requested state.
    /// @param state convert to this state.
    void
    convert_to_state(State state);

    /// Return the reference particle
    Reference_particle &
    get_reference_particle();

    Reference_particle const&
    get_reference_particle() const;

    /// Get the array containing the macroparticles on this processor.
    /// The array has length (length,7), where length of the array may be
    /// larger local_num. The macroparticle state vectors are stored in
    /// array[0:local_num,0:6] and the macroparticle IDs are stored in
    /// array[0:local_num,6]. Use get_local_num() to obtain local_num.
    MArray2d_ref
    get_local_particles();

    Const_MArray2d_ref
    get_local_particles() const;

    /// Get the particle charge in units of e.
    int
    get_particle_charge() const;

    /// Get the particle mass in units GeV/c^2.
    double
    get_mass() const;

    /// Get the real number of particles represented by the bunch.
    double
    get_real_num() const;

    /// Get the number of macroparticles stored on this processor.
    int
    get_local_num() const;

    /// Get the total number of macroparticles.
    int
    get_total_num() const;

    /// Get the (fixed-t or fixed-z) state.
    State
    get_state() const;

    /// Get the communicator
    Commxx const&
    get_comm() const;

    /// Add a copy of the particles in bunch to the current bunch. The
    /// injected bunch must have the same macroparticle weight, i.e.,
    /// real_num/total_num. If the state vectors of the reference particles
    /// of the two bunches differ, the particles will be shifted accordingsly.
    void
    inject(Bunch const& bunch);

    virtual
    ~Bunch();
};

#endif /* BUNCH_H_ */
