#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <string>
#include <boost/shared_ptr.hpp>
#include "synergia/bunch/bunch.h"
//#include "synergia/foundation/generalized_diagnostics.h"
#include "synergia/foundation/diagnostics_write_helper.h"
#include "synergia/utils/hdf5_serial_writer.h"

/// Diagnostics is an abstract base class for bunch diagnostics classes
struct Core_diagnostics
{
    static MArray1d
    calculate_mean(Bunch const& bunch);

    static double
    calculate_z_mean(Bunch const& bunch);

    static MArray1d
    calculate_std(Bunch const& bunch, MArray1d_ref const& mean);

    static MArray2d
    calculate_mom2(Bunch const& bunch, MArray1d_ref const& mean);

    static MArray1d
    calculate_min(Bunch const& bunch);

    static MArray1d
    calculate_max(Bunch const& bunch);
};

class Diagnostics
{
private:
    std::string name;
    std::string filename;
    Bunch_sptr bunch_sptr;
    bool have_bunch_;
    Diagnostics_write_helper * write_helper_ptr;
    bool have_write_helper_;

public:

    Diagnostics(std::string const& name, std::string const& filename);

    // Default constructor for serialization use only
    Diagnostics();

    virtual std::string const&
    get_filename() const;

    virtual void
    set_bunch_sptr(Bunch_sptr bunch_sptr);

    virtual bool
    have_bunch() const;

    virtual void
    delete_write_helper_ptr();

    virtual Diagnostics_write_helper *
    new_write_helper_ptr();

    virtual bool
    have_write_helper() const;

    virtual Diagnostics_write_helper *
    get_write_helper_ptr();

    Bunch &
    get_bunch();
    /// Multiple serial diagnostics can be written to a single file.
    virtual bool
    is_serial() const = 0;
    /// Update the diagnostics
    virtual void
    update() = 0;
    /// Write the diagnostics to the file
    virtual void
    write() = 0;
    /// Update the diagnostics and write them to the file
    virtual void
    update_and_write()
    {
        update();
        write();
    }
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(name);
            ar & BOOST_SERIALIZATION_NVP(filename);
            ar & BOOST_SERIALIZATION_NVP(bunch_sptr);
            ar & BOOST_SERIALIZATION_NVP(have_bunch_);
            ar & BOOST_SERIALIZATION_NVP(write_helper_ptr);
            ar & BOOST_SERIALIZATION_NVP(have_write_helper_);
        }
    virtual
    ~Diagnostics();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics)
typedef boost::shared_ptr<Diagnostics > Diagnostics_sptr;

/// Diagnostics_basic provides the minimal set of statistical
/// quantities to be calculated for a Bunch.
class Diagnostics_basic : public Diagnostics
{
private:
    bool have_writers;
    double s;
    Hdf5_serial_writer<double > * writer_s;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double trajectory_length;
    Hdf5_serial_writer<double > * writer_trajectory_length;
    int num_particles;
    Hdf5_serial_writer<int > * writer_num_particles;
    double real_num_particles;
    Hdf5_serial_writer<double > * writer_real_num_particles;
    MArray1d mean;
    Hdf5_serial_writer<MArray1d_ref > * writer_mean;
    MArray1d std;
    Hdf5_serial_writer<MArray1d_ref > * writer_std;
    MArray1d min;
    Hdf5_serial_writer<MArray1d_ref > * writer_min;
    MArray1d max;
    Hdf5_serial_writer<MArray1d_ref > * writer_max;


public:
    /// Create a Diagnostics_basic object
    /// @param filename filename for output
    Diagnostics_basic(std::string const& filename);

    // Default constructor for serialization use only
    Diagnostics_basic();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_basic class is serial.
    virtual bool
    is_serial() const;

    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

    /// Update the diagnostics
    virtual void
    update();

    /// Get the distance from the origin along the reference trajectory in
    /// meters.
    virtual double
    get_s() const;

    /// Get the number of complete repetitions.
    virtual int
    get_repetition() const;

    /// Get the total distance along the reference trajectory in meters.
    virtual double
    get_trajectory_length() const;

    /// Get the total number of macroparticles in the bunch
    virtual int
    get_num_particles() const;

    /// Get the total number of real particles represented by the bunch
    virtual double
    get_real_num_particles() const;

    /// Get a six-dimensional vector of the means of each phase-space
    /// coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref
    get_mean() const;

    /// Get a six-dimensional vector of the standard deviations of each
    /// phase-space coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref
    get_std() const;

    virtual const MArray1d
    get_min() const;

    virtual const MArray1d
    get_max() const;

    virtual void
    write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
            ar & BOOST_SERIALIZATION_NVP(have_writers);
            ar & BOOST_SERIALIZATION_NVP(s);
            ar & BOOST_SERIALIZATION_NVP(writer_s);
            ar & BOOST_SERIALIZATION_NVP(repetition);
            ar & BOOST_SERIALIZATION_NVP(writer_repetition);
            ar & BOOST_SERIALIZATION_NVP(trajectory_length);
            ar & BOOST_SERIALIZATION_NVP(writer_trajectory_length);
            ar & BOOST_SERIALIZATION_NVP(num_particles);
            ar & BOOST_SERIALIZATION_NVP(writer_num_particles);
            ar & BOOST_SERIALIZATION_NVP(real_num_particles);
            ar & BOOST_SERIALIZATION_NVP(writer_real_num_particles);
            ar & BOOST_SERIALIZATION_NVP(mean);
            ar & BOOST_SERIALIZATION_NVP(writer_mean);
            ar & BOOST_SERIALIZATION_NVP(std);
            ar & BOOST_SERIALIZATION_NVP(writer_std);
            ar & BOOST_SERIALIZATION_NVP(min);
            ar & BOOST_SERIALIZATION_NVP(writer_min);
            ar & BOOST_SERIALIZATION_NVP(max);
            ar & BOOST_SERIALIZATION_NVP(writer_max);
        }

    virtual
    ~Diagnostics_basic();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_basic)
typedef boost::shared_ptr<Diagnostics_basic > Diagnostics_basic_sptr;

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_full2 : public Diagnostics
{
private:
    bool have_writers;
    double s;
    Hdf5_serial_writer<double > * writer_s;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double trajectory_length;
    Hdf5_serial_writer<double > * writer_trajectory_length;
    int num_particles;
    Hdf5_serial_writer<int > * writer_num_particles;
    double real_num_particles;
    Hdf5_serial_writer<double > * writer_real_num_particles;
    MArray1d mean;
    Hdf5_serial_writer<MArray1d_ref > * writer_mean;
    MArray1d std;
    Hdf5_serial_writer<MArray1d_ref > * writer_std;
    MArray1d min;
    Hdf5_serial_writer<MArray1d_ref > * writer_min;
    MArray1d max;
    Hdf5_serial_writer<MArray1d_ref > * writer_max;
    MArray2d mom2;
    Hdf5_serial_writer<MArray2d_ref > * writer_mom2;
    MArray2d corr;
    Hdf5_serial_writer<MArray2d_ref > * writer_corr;
    double emitx, emity, emitz, emitxy, emitxyz;
    Hdf5_serial_writer<double > *writer_emitx, *writer_emity, *writer_emitz,
            *writer_emitxy, *writer_emitxyz;
    virtual void
    update_full2();
    virtual void
    update_emittances();

public:
    /// Create a Diagnostics_basic object
    /// @param bunch the Bunch
    /// @param filename filename for output
    Diagnostics_full2(std::string const& filename);

    // Default constructor for serialization use only
    Diagnostics_full2();

    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_full2 class is serial.
    virtual bool
    is_serial() const;

    /// Update the diagnostics
    virtual void
    update();

    /// Get the distance from the origin along the reference trajectory in
    /// meters.
    virtual double
    get_s() const;

    /// Get the number of complete repetitions.
    virtual int
    get_repetition() const;

    /// Get the total distance along the reference trajectory in meters.
    virtual double
    get_trajectory_length() const;

    /// Get the total number of macroparticles in the bunch
    virtual int
    get_num_particles() const;

    /// Get the total number of real particles represented by the bunch
    virtual double
    get_real_num_particles() const;

    /// Get a six-dimensional vector of the means of each phase-space
    /// coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref
    get_mean() const;

    /// Get a six-dimensional vector of the standard deviations of each
    /// phase-space coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref
    get_std() const;

    virtual const MArray1d
    get_min() const;

    virtual const MArray1d
    get_max() const;

    /// Get a 6x6 matrix of the second moments of the phase-space coordinates.
    /// The units are Synergia units.
    virtual Const_MArray2d_ref
    get_mom2() const;

    /// Get a 6x6 matrix of the correlation coefficients of the phase-space
    /// coordinates.
    virtual Const_MArray2d_ref
    get_corr() const;

    /// Get the horizontal emittance.
    /// Currently reported in unnatural Synergia units.
    virtual double
    get_emitx() const;

    /// Get the vertical emittance.
    /// Currently reported in unnatural Synergia units.
    virtual double
    get_emity() const;

    /// Get the longitudinal emittance.
    /// Currently reported in unnatural Synergia units.
    virtual double
    get_emitz() const;

    /// Get the (4D) transverse emittance.
    /// Currently reported in unnatural Synergia units.
    virtual double
    get_emitxy() const;

    /// Get the (6D) full emittance.
    /// Currently reported in unnatural Synergia units.
    virtual double
    get_emitxyz() const;

    virtual void
    write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                    & BOOST_SERIALIZATION_NVP(have_writers)
                    & BOOST_SERIALIZATION_NVP(s)
                    & BOOST_SERIALIZATION_NVP(writer_s)
                    & BOOST_SERIALIZATION_NVP(repetition)
                    & BOOST_SERIALIZATION_NVP(writer_repetition)
                    & BOOST_SERIALIZATION_NVP(trajectory_length)
                    & BOOST_SERIALIZATION_NVP(writer_trajectory_length)
                    & BOOST_SERIALIZATION_NVP(num_particles)
                    & BOOST_SERIALIZATION_NVP(writer_num_particles)
                    & BOOST_SERIALIZATION_NVP(real_num_particles)
                    & BOOST_SERIALIZATION_NVP(writer_real_num_particles)
                    & BOOST_SERIALIZATION_NVP(mean)
                    & BOOST_SERIALIZATION_NVP(writer_mean)
                    & BOOST_SERIALIZATION_NVP(std)
                    & BOOST_SERIALIZATION_NVP(writer_std)
                    & BOOST_SERIALIZATION_NVP(min)
                    & BOOST_SERIALIZATION_NVP(writer_min)
                    & BOOST_SERIALIZATION_NVP(max)
                    & BOOST_SERIALIZATION_NVP(writer_max)
                    & BOOST_SERIALIZATION_NVP(mom2)
                    & BOOST_SERIALIZATION_NVP(writer_mom2)
                    & BOOST_SERIALIZATION_NVP(corr)
                    & BOOST_SERIALIZATION_NVP(writer_corr)
                    & BOOST_SERIALIZATION_NVP(emitx)
                    & BOOST_SERIALIZATION_NVP(emity)
                    & BOOST_SERIALIZATION_NVP(emitz)
                    & BOOST_SERIALIZATION_NVP(emitxy)
                    & BOOST_SERIALIZATION_NVP(emitxyz)
                    & BOOST_SERIALIZATION_NVP(writer_emitx)
                    & BOOST_SERIALIZATION_NVP(writer_emity)
                    & BOOST_SERIALIZATION_NVP(writer_emitz)
                    & BOOST_SERIALIZATION_NVP(writer_emitxy)
                    & BOOST_SERIALIZATION_NVP(writer_emitxyz);
        }

    virtual
    ~Diagnostics_full2();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_full2)
typedef boost::shared_ptr<Diagnostics_full2 > Diagnostics_full2_sptr;

/// Diagnostics_particles dumps the state of particles in a bunch
class Diagnostics_particles : public Diagnostics
{
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
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
            ar & BOOST_SERIALIZATION_NVP(have_writers);
            ar & BOOST_SERIALIZATION_NVP(max_particle_id);
        }

    virtual
    ~Diagnostics_particles();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_particles)
typedef boost::shared_ptr<Diagnostics_particles > Diagnostics_particles_sptr;

/// Diagnostics_track records the phase space coordinates of a single particle.
/// Particles will only be tracked if they stay on the same processor.
/// Lost particles that are somehow restored or particles not available when
/// the first update is called will also not be tracked.
class Diagnostics_track : public Diagnostics
{
private:
    bool have_writers;
    bool found;
    bool first_search;
    int last_index;
    int particle_id;
    double s;
    Hdf5_serial_writer<double > * writer_s;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double trajectory_length;
    Hdf5_serial_writer<double > * writer_trajectory_length;
    MArray1d coords;
    Hdf5_serial_writer<MArray1d_ref > * writer_coords;
    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

public:
    /// Create an empty Diagnostics_track object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param particle_id the particle ID to track
    Diagnostics_track(std::string const& filename,
            int particle_id);

    // Default constructor for serialization use only
    Diagnostics_track();

    virtual Diagnostics_write_helper *
    new_write_helper_ptr();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_full2 class is serial.
    virtual bool
    is_serial() const;

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                    & BOOST_SERIALIZATION_NVP(have_writers)
                    & BOOST_SERIALIZATION_NVP(found)
                    & BOOST_SERIALIZATION_NVP(first_search)
                    & BOOST_SERIALIZATION_NVP(last_index)
                    & BOOST_SERIALIZATION_NVP(particle_id)
                    & BOOST_SERIALIZATION_NVP(s)
                    & BOOST_SERIALIZATION_NVP(writer_s)
                    & BOOST_SERIALIZATION_NVP(repetition)
                    & BOOST_SERIALIZATION_NVP(writer_repetition)
                    & BOOST_SERIALIZATION_NVP(trajectory_length)
                    & BOOST_SERIALIZATION_NVP(writer_trajectory_length)
                    & BOOST_SERIALIZATION_NVP(coords)
                    & BOOST_SERIALIZATION_NVP(writer_coords);
        }

    virtual
    ~Diagnostics_track();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_track)
typedef boost::shared_ptr<Diagnostics_track > Diagnostics_track_sptr;

class Diagnostics_reference_particle : public Diagnostics
{
private:
    bool have_writers;
    Hdf5_serial_writer<double > * writer_beta;
    Hdf5_serial_writer<double > * writer_gamma;
    Hdf5_serial_writer<MArray1d_ref > * writer_state;
    Hdf5_serial_writer<double > * writer_s;

public:
    /// Create a Diagnostics_reference_particle object
    /// @param bunch the Bunch
    /// @param filename filename for output
    Diagnostics_reference_particle(
            std::string const& filename);

    // Default constructor for serialization use only
    Diagnostics_reference_particle();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_reference_particle class is serial.
    virtual bool
    is_serial() const;

    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                    & BOOST_SERIALIZATION_NVP(have_writers)
                    & BOOST_SERIALIZATION_NVP(writer_beta)
                    & BOOST_SERIALIZATION_NVP(writer_gamma)
                    & BOOST_SERIALIZATION_NVP(writer_state)
                    & BOOST_SERIALIZATION_NVP(writer_s);
        }

    virtual
    ~Diagnostics_reference_particle();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_reference_particle)
typedef boost::shared_ptr<Diagnostics_reference_particle >
        Diagnostics_reference_particle_sptr;

#endif /* DIAGNOSTICS_H_ */
