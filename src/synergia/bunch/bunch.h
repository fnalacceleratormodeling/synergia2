#ifndef BUNCH_H_
#define BUNCH_H_

#include <sstream>
#include <vector>
#include <memory>

#include "synergia/foundation/reference_particle.h"

#include "synergia/bunch/bunch_particles.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/diagnostics_loss.h"

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/restrict_extension.h"
#include "synergia/utils/logger.h"

enum class ParticleGroup
{
    regular = 0,
    spectator = 1
};

enum class LongitudinalBoundary
{
    open = 0,
    periodic = 1,
    aperture = 2,
    bucket_barrier = 3
};

/// Represents a macroparticle bunch distributed across the processors
/// in a comm_sptrunicator.
class Bunch
{
private:

    using PG = ParticleGroup;
    using LB = LongitudinalBoundary;

public:
    /*! \enum State The state of the bunch is captured at a fixed  s (or z, longitudinal coordinate)
     or at a fixed time.  In the former case, particles are found within a range of different time
     coordinates while in the later case particles position along the beam axis do vary.
     A change of state is accomplish via the fixed_t_z_converter class.
     */
    constexpr static const int x    = 0;
    constexpr static const int xp   = 1;
    constexpr static const int y    = 2;
    constexpr static const int yp   = 3;
    constexpr static const int z    = 4;
    constexpr static const int zp   = 5;
    constexpr static const int cdt  = 4;
    constexpr static const int dpop = 5;
    constexpr static const int id   = 6;

    constexpr static const int particle_index_null = 
        BunchParticles::particle_index_null;

private:

    Commxx comm;    

    int bucket_index;

    // meaning of bounary_param for each boundary condition:
    // open (N/A), periodic (z-period), z-cut (longitudinal_extent),
    // bucket_barrier (bucket_length)
    LongitudinalBoundary boundary;
    double boundary_param;

    // reference particle and design reference particle
    Reference_particle ref_part;
    Reference_particle design_ref_part;

    int particle_charge;

    double real_num;

    // parts[0]: PG::regular particles
    // parts[1]: spectator particles
    std::array<BunchParticles, 2> parts;

    // diagnostics
    std::map<std::string, std::unique_ptr<Diagnostics>> diags;

    // diagnostics for particle losses
    std::unique_ptr<Diagnostics_loss> diag_aperture;
    std::unique_ptr<Diagnostics_loss> diag_zcut;


private:

    std::string get_local_particles_serialization_path() const;

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
    Bunch(  Reference_particle const& reference_particle, 
            int total_num,
            double real_num, 
            Commxx comm,
            int total_spectator_num = 0 );
    
    // non-copyable but moveable
    Bunch(Bunch const&) = delete;
    Bunch(Bunch &&) = default;

    ///
    /// Set the particle charge
    /// @param particle_charge in units of e.
    void set_particle_charge(int charge)
    { particle_charge = charge; }

    ///
    /// Set the number of real particles represented by the bunch.
    /// @param real_num the new real number of particles
    void set_real_num(double num)
    { real_num = num; }

    /// Return the reference particle
    Reference_particle &      get_reference_particle()       
    { return ref_part; }

    Reference_particle const& get_reference_particle() const 
    { return ref_part; }

    Reference_particle &      get_design_reference_particle()       
    { return design_ref_part; }

    Reference_particle const& get_design_reference_particle() const 
    { return design_ref_part; }

    void set_design_reference_particle(Reference_particle const & ref_part)
    { design_ref_part = ref_part; }

    // BunchParticles
    BunchParticles &      get_bunch_particles(ParticleGroup pg = PG::regular)
    { return parts[(int)pg]; }

    BunchParticles const& get_bunch_particles(ParticleGroup pg = PG::regular) const
    { return parts[(int)pg]; }

    /// Get the array containing the macroparticles on this processor.
    /// The array has length (length,7), where length of the array may be
    /// larger local_num. The macroparticle state vectors are stored in
    /// array[0:local_num,0:6] and the macroparticle IDs are stored in
    /// array[0:local_num,6]. Use get_local_num() to obtain local_num.
    Particles          get_local_particles(ParticleGroup pg = PG::regular)       
    { return get_bunch_particles(pg).parts; }

    ConstParticles     get_local_particles(ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).parts; }

    HostParticles      get_host_particles(ParticleGroup pg = PG::regular)       
    { return get_bunch_particles(pg).hparts; }

    ConstHostParticles get_host_particles(ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).hparts; }

    /// getters of particle array dimensions
    int get_total_num(ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).total_num(); }

    int get_local_num(ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).local_num(); }

    int get_local_num_aligned (ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).local_num_aligned(); }

    int get_local_num_padded  (ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).local_num_padded(); }

    int get_local_num_slots   (ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).local_num_slots(); }

    int get_local_num_padding (ParticleGroup pg = PG::regular) const 
    { return get_bunch_particles(pg).local_num_padding(); }

    int get_local_num_lost    (ParticleGroup pg = PG::regular) const 
    { return get_local_num_slots(pg) - get_local_num_padded(pg); }

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
    void set_local_num(int num, ParticleGroup pg = PG::regular)
    { get_bunch_particles(pg).set_local_num(num); }

    ///
    /// Set the total number (and the real number) of particles
    void set_total_num(int num, ParticleGroup pg = PG::regular)
    { get_bunch_particles(pg).set_total_num(num); }
 
    ///
    /// Update the total number and real number of particles after the local
    /// number has been changed. Requires comm_sptrunication.
    void update_total_num()
    { 
        get_bunch_particles(PG::spectator).update_total_num(); 

        auto & bp = get_bunch_particles(PG::regular);
        int old_total = bp.update_total_num();
        real_num = old_total ? bp.total_num() * real_num / old_total : 0.0;
    }

    // aperture operation
    // discard the particles (by moving them to the tail of the array) filtered out
    // by the aperture, returns the number of particles discarded
    template<typename AP>
    int apply_aperture(AP const& ap, ParticleGroup pg = PG::regular);
  
    // checkout (deep_copy) the entire particle array from device
    // memory to the host memory for user to access the latest
    // particle data
    void checkout_particles(ParticleGroup pg = PG::regular) 
    { get_bunch_particles(pg).checkout_particles(); }

    void checkin_particles (ParticleGroup pg = PG::regular)
    { get_bunch_particles(pg).checkin_particles(); }

    // checkout (deep_copy) num particles starting from idx, and
    // store them in a host array
    // TODO: when compiled for host, it returns a subview to the original
    // particle data -- overhead for the operation should be minimal
    karray2d_row 
    get_particles_in_range(int idx, int num, ParticleGroup pg = PG::regular) const
    { return get_bunch_particles(pg).get_particles_in_range(idx, num); }

    karray1d_row 
    get_particle(int idx, ParticleGroup pg = PG::regular) const
    { return get_bunch_particles(pg).get_particle(idx); }

    // find the index of the given particle_id (pid)
    // if last_idx is provided, it does the search form the last_idx first
    // returns particle_index_null if the given particle_id is not found
    int search_particle(int pid, int last_idx = particle_index_null, 
            ParticleGroup pg = PG::regular) const
    { return get_bunch_particles(pg).search_particle(pid, last_idx); }

    void print_particle(size_t idx, Logger& logger, 
            ParticleGroup pg = PG::regular) const
    { get_bunch_particles(pg).print_particle(idx, logger); }

#if 0
    std::array<size_t, 2> get_particle_strides(ParticleGroup pg = PG::regular) const
    { 
        std::array<size_t, 2> strides;
        get_bunch_particles(pg).parts.stride(strides.data());
        return strides;
    }
#endif

    /// Get the particle charge in units of e.
    int get_particle_charge() const { return particle_charge; }

    /// Get the particle mass in units GeV/c^2.
    double get_mass() const { return ref_part.get_four_momentum().get_mass(); }

    /// Get the real number of particles represented by the bunch.
    double get_real_num() const { return real_num; }

    /// longitudinal boundary conditions
    void set_longitudinal_bounadry(LB lb, double param = 0.0)
    { boundary = lb; boundary_param = param; }

    std::pair<LB, double> get_longitudinal_bounadry() const
    { return std::make_pair(boundary, boundary_param); }

    // bucket index
    void set_bucket_index(int index)      { bucket_index = index; }
    int  get_bucket_index() const         { return bucket_index; }
    bool is_bucket_index_assigned() const { return bucket_index >= 0; }

    /// Get the communicator
    Commxx const& get_comm() const { return comm; }

    // Diagnostics
    template<typename Diag>
    void add_diagnostics(std::string const& name, Diag&& diag)
    { diags.emplace(name, std::make_unique<Diag>(std::move(diag)))
        .first->second->set_bunch(*this); }

    Diagnostics & get_diag(std::string const & name);

    void set_diag_loss_aperture(Diagnostics_loss && diag)
    { diag_aperture = std::make_unique<Diagnostics_loss>(std::move(diag)); }

    void set_diag_loss_zcut(Diagnostics_loss && diag)
    { diag_zcut = std::make_unique<Diagnostics_loss>(std::move(diag)); }

    /// Add a copy of the particles in bunch to the current bunch. The
    /// injected bunch must have the same macroparticle weight, i.e.,
    /// real_num/total_num. If the state vectors of the reference particles
    /// of the two bunches differ, the particles will be shifted accordingly.
    void inject(Bunch const& bunch);
    
    void read_file(std::string const &);

    void check_pz2_positive()
    {
        get_bunch_particles(ParticleGroup::regular).check_pz2_positive();
        get_bunch_particles(ParticleGroup::spectator).check_pz2_positive();
    }
    
    /// serialization
    template<class Archive> 
    void save(Archive & ar, const unsigned int version) const;

    template<class Archive>
    void load(Archive & ar, const unsigned int version);
};

template<typename AP>
inline int Bunch::apply_aperture(AP const& ap, ParticleGroup pg)
{ 
    int ndiscarded = get_bunch_particles(pg).apply_aperture(ap); 

    if (ndiscarded && diag_aperture)
    {
        auto const& bp = get_bunch_particles(pg);
        auto discarded = get_particles_in_range(bp.local_num_padded(), bp.last_discarded(), pg);
        //diag_aperture->update(discarded);
        //diag_aperture->write();
    }

    return ndiscarded;
}


#endif /* BUNCH_H_ */
