#ifndef BUNCH_H_
#define BUNCH_H_

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/commxx.h"
#include "synergia/bunch/fixed_t_z_converter.h"
#include "boost/shared_ptr.hpp"
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/restrict_extension.h"

#include <mpi.h>

#include <string>
#include <utility>
#include <vector>

/// Represents a macroparticle bunch distributed across the processors
/// in a comm_sptrunicator.
class Bunch
{
public:
    /*! \enum State The state of the bunch is captured at a fixed  s (or z, longitudinal coordinate)
     or at a fixed time.  In the former case, particles are found within a range of different time
     coordinates while in the later case particles position along the beam axis do vary.
     A change of state is accomplish via the fixed_t_z_converter class.
     */
    enum State
    {   
        fixed_z_lab = 1,
        fixed_t_lab = 2,
        fixed_t_bunch = 3,       
        fixed_z_bunch = 4
    };

    static const int x = 0;
    static const int xp = 1;
    static const int y = 2;
    static const int yp = 3;
    static const int z = 4;
    static const int zp = 5;
    static const int cdt = 4;
    static const int dpop = 5;
    static const int id = 6;

    const static int particle_alignment;

    // longitudinal boundary conditions
    enum LongitudinalBoundary
    {
        lb_open = 0,
        lb_periodic = 1,
        lb_aperture = 2,
        lb_bucket_barrier = 3
    };

private:

    LongitudinalBoundary boundary = lb_open;
    double boundary_param = 0.0;  // NA, z-period, longitudinal_extent, or bucket_length

    Reference_particle reference_particle;
    Reference_particle design_reference_particle;
    int particle_charge = 1;

    /*
     * Local Particle Array Memory Layout:
     *
     *   P: particle, O: padding, L: lost particle
     *
     *   +=====+
     *   |  P  |
     *   +-----+
     *   |  P  |
     *   +-----+
     *   | ... |
     *   +=====+  <- local_num
     *   |  O  |
     *   +-----+  <- local_num_aligned
     *   |  O  |
     *   +-----+
     *   | ... |
     *   +=====+  <- local_num_padded
     *   |  L  |
     *   +-----+
     *   |  L  |
     *   +-----+
     *   | ... |
     *   +=====+  <- local_num_slots
     *
     *   * number of padding slots  = local_num_padded - local_num
     *   * number of lost particles = local_num_slots - local_num_padded
     *
     *   At bunch construction the size of padding (num_padding) is decided 
     *   such that the local_num_slots is always aligned (depending on the 
     *   vector specification, e.g., SSE or AVX or AVX512). 
     *
     *   local_num_aligned is initialized in the range [local_num, 
     *   local_num_paded], and gets adjusted everytime the local_num 
     *   changes such tht local_num_aligned is always aligned.
     *
     */

    int local_num = 0;
    int local_num_aligned, local_num_padded, local_num_slots;
    int local_s_num, local_s_num_aligned, local_s_num_padded, local_s_num_slots;

    int total_num, total_s_num;

    double real_num;

    double * storage = nullptr;
    double * s_storage = nullptr;

    MArray2d_ref *local_particles = nullptr;
    MArray2d_ref *local_s_particles = nullptr;

    int bucket_index;
    bool bucket_index_assigned;

    int sort_period, sort_counter;

    State state;
    Commxx_sptr comm_sptr;    
    Fixed_t_z_zeroth default_converter;
    Fixed_t_z_converter *converter_ptr;
    // Fixed_t_z_alex default_converter;
    //  Fixed_t_z_synergia20 default_converter;
    void
    assign_ids(int local_offset);

    void
    assign_spectator_ids(int local_offset);

    std::string
    get_local_particles_serialization_path() const;

    void
    construct(int total_num, double real_num, int total_s_num);

public:
    //!
    //! Constructor:
    //! Allocates memory for the particles and assigns particle ID's,
    //!    but does not fill the phase space values in any way.
    //!
    //! To fill the bunch with particles, use the populate methods.
    /// @param reference_particle the reference particle for the bunch.
    /// @param total_num the total number of macroparticles in the bunch
    /// @param real_num the number of real particles represented by the bunch.
    /// @param bucket_index the bucket number the  bunch occupies, used for multi-bunch simulations
    /// @param comm_sptr the comm_sptrunicator.
    Bunch(Reference_particle const& reference_particle, int total_num,
            double real_num, Commxx_sptr comm_sptr);

    Bunch(Reference_particle const& reference_particle, int total_num, int total_spectator_num,
            double real_num, Commxx_sptr comm_sptr);

     ///// Obsolete, please replace the following constructor with the previous one followed by 
     /////set_particle_charge(particle_charge)    
     //    Bunch(Reference_particle const& reference_particle, int total_num,
     //            double real_num, Commxx_sptr comm_sptr, int particle_charge);

    /////Obsolete, please replace the following constructor with the previous one followed by 
    /////set_z_period_length(z_period_length)
    /////set_bucket_index(bucket_index)
    //Bunch(Reference_particle const& reference_particle, int total_num,
    //        double real_num, Commxx_sptr comm_sptr, double z_period_length,
    //        int bucket_index = 0);

    
    
    /// Default constructor for serialization use only
    Bunch();

    //!
    //! Copy constructor
    Bunch(Bunch const& bunch);
    //!
    //! Assignment constructor
    Bunch &
    operator=(Bunch const& bunch);

    ~Bunch();

    ///
    /// Set the particle charge
    /// @param particle_charge in units of e.
    void
    set_particle_charge(int particle_charge);

    ///
    /// Set the number of real particles represented by the bunch.
    /// @param real_num the new real number of particles
    void
    set_real_num(double real_num);

    ///
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

    void
    set_local_spectator_num(int local_s_num);

    void
    expand_local_num(int num, int added_lost);

    void
    expand_local_spectator_num(int s_num, int added_lost);

    ///
    /// Update the total number and real number of particles after the local
    /// number has been changed. Requires comm_sptrunication.
    void
    update_total_num();
    
    ///
    /// Set the total number (and the real number) of particles
    void
    set_total_num(int totalnum);
    
    /// For a given number of particles, returns the next alignement position
    /// @param num number of particles
    int
    calculate_aligned_pos(int num);

    /// For a given number of particles, returns the needed size of padding
    /// for the particle array to be aligned
    /// @param num number of particles
    int
    calculate_padding_size(int num);

    ///
    /// Set the period for periodic_sort and reset the counter
    /// Periods less than zero will prohibit sorting.
    /// @param period
    void
    set_sort_period(int period);

    /// Sort the particles
    /// @param index the particle index on which to sort
    void
    sort(int index);

    /// Sort the particles every sort period calls
    /// @param index the particle index on which to sort
    void
    periodic_sort(int index);

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

    Reference_particle &
    get_design_reference_particle();

    Reference_particle const &
    get_design_reference_particle() const;

    void
    set_design_reference_particle(Reference_particle const & ref_part);

    /// Get the array containing the macroparticles on this processor.
    /// The array has length (length,7), where length of the array may be
    /// larger local_num. The macroparticle state vectors are stored in
    /// array[0:local_num,0:6] and the macroparticle IDs are stored in
    /// array[0:local_num,6]. Use get_local_num() to obtain local_num.
    MArray2d_ref
    get_local_particles();

    Const_MArray2d_ref
    get_local_particles() const;

    MArray2d_ref
    get_local_spectator_particles();

    Const_MArray2d_ref
    get_local_spectator_particles() const;

    /// Get the particle charge in units of e.
    int
    get_particle_charge() const;

    /// Get the particle mass in units GeV/c^2.
    double
    get_mass() const;

    /// Get the real number of particles represented by the bunch.
    double
    get_real_num() const;

    /// Get the period length of the bunch
    double
    get_z_period_length() const;
    
    /// Set the period length of the bunch and make the bunch z_periodic     
    void
    set_z_period_length(double z_period_length) ;
       
    /// Is the bunch periodic?
    bool
    is_z_periodic() const; 
    
    /// Get the the bunch extent if the longitudinal aperture is present
    double
    get_longitudinal_aperture_length() const;

    /// Set the longitudinal_extent of the bunch and make the longitudinal aperture true
    void
    set_longitudinal_aperture_length(double longitudinal_length);
    
    /// True when the longitudinal aperture is present and the bunch is cut after every  operation
    /// longitudinally outside the extent [-longitudinal_extent/2,  longitudinal_extent/2]    
    bool
    has_longitudinal_aperture() const; 

    double
    get_bucket_barrier_length() const;

    void
    set_bucket_barrier_length(double length);

    bool
    has_bucket_barrier() const;

    void
    set_open_longitudinal_boundary();

    bool
    is_open_longitudinal_boundary() const;

    /// A general interface to set the bunch boundary condition
    void
    set_longitudinal_boundary(LongitudinalBoundary boundary, double period = 0.0);

    /// The general interface for querying the bunch boundary condition
    std::pair<Bunch::LongitudinalBoundary, double>
    get_longitudinal_boundary() const;

    /// Get the number of macroparticles stored on this processor.
    int
    get_local_num() const;

    /// Get the number of padded macroparticles (first dimension of the particles[][] array)
    int
    get_local_num_slots() const;

    int
    get_local_num_aligned() const;

    int
    get_local_num_padded() const;

    int
    get_local_num_padding() const;

    int
    get_local_num_lost() const;

    /// Get the total number of macroparticles.
    int
    get_total_num() const;

    /// Get the number of spectator particles stored on this processor.
    int
    get_local_spectator_num() const;

    int
    get_local_spectator_num_slots() const;

    int
    get_local_spectator_num_aligned() const;

    int
    get_local_spectator_num_padded() const;

    int
    get_local_spectator_num_padding() const;

    int
    get_local_spectator_num_lost() const;

    /// Get the total number of spectator particles.
    int
    get_total_spectator_num() const;

    /// Get the period for periodic_sort
    int
    get_sort_period() const;

    void
    set_bucket_index(int index);

    std::size_t
    get_bucket_index() const;
    
    bool
    is_bucket_index_assigned() const;

    /// Get the (fixed-t or fixed-z) state.
    State
    get_state() const;

    /// Get the communicator
    Commxx const&
    get_comm() const;

    /// Get the communicator
    Commxx_sptr
    get_comm_sptr() const;

    /// Add a copy of the particles in bunch to the current bunch. The
    /// injected bunch must have the same macroparticle weight, i.e.,
    /// real_num/total_num. If the state vectors of the reference particles
    /// of the two bunches differ, the particles will be shifted accordingly.
    void
    inject(Bunch const& bunch);
    
    void
    read_file(std::string const &);

    void
    check_pz2_positive();
    
    void set_arrays(double * RESTRICT &xa, double * RESTRICT &xpa,
                    double * RESTRICT &ya, double * RESTRICT &ypa,
                    double * RESTRICT &cdta, double * RESTRICT &dpopa);

    void set_spectator_arrays(double * RESTRICT &xa, double * RESTRICT &xpa,
                    double * RESTRICT &ya, double * RESTRICT &ypa,
                    double * RESTRICT &cdta, double * RESTRICT &dpopa);

    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const;

    template<class Archive>
        void
        load(Archive & ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

using Bunch_sptr = boost::shared_ptr<Bunch>;
using Bunches = std::vector<Bunch_sptr>;

#endif /* BUNCH_H_ */
