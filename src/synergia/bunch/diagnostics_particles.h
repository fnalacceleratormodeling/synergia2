#ifndef DIAGNOSTICS_PARTICLES_H_
#define DIAGNOSTICS_PARTICLES_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_particles dumps the state of particles in a bunch
class Diagnostics_particles : public Diagnostics
{
public:
    static const char name[];
private:
    bool have_writers;
    int min_particle_id, max_particle_id;
    void
    receive_other_local_particles(std::vector<int > const& local_nums,
            Hdf5_file_sptr file_sptr);
    void
    send_local_particles();
public:
    /// Create a Diagnostics_particles object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param min_particle_id the lowest particle id to write (defaults to 0)
    /// @param max_particle_id the highest particle id to write (0 indicates no limit, hence min,max = 0,0 writes all particles)
    /// @param write_skip write every write_skip turns
    Diagnostics_particles(std::string const& filename, int min_particle_id = 0,
            int max_particle_id = 0, int write_skip = 1);

    // Default constructor for serialization use only
    Diagnostics_particles();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_particles class is not serial.
    virtual bool
    is_serial() const;

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Diagnostics_particles();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_particles)
typedef boost::shared_ptr<Diagnostics_particles > Diagnostics_particles_sptr;

#endif /* DIAGNOSTICS_PARTICLES_H_ */
