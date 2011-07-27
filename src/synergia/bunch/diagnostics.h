#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <string>
#include <boost/shared_ptr.hpp>

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics_write_helper.h"
#include "synergia/utils/hdf5_serial_writer.h"

/// Diagnostics is an abstract base class for Diagnostics classes
class Diagnostics
{

public:
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
    update_and_write();

    static MArray1d
    calculate_mean(Bunch const& bunch);

    static MArray1d
    calculate_std(Bunch const& bunch, MArray1d_ref const& mean);

    static MArray2d
    calculate_mom2(Bunch const& bunch, MArray1d_ref const& mean);

    static MArray1d
    calculate_bunchmin(Bunch const& bunch);

    static MArray1d
    calculate_bunchmax(Bunch const& bunch);

    virtual
    ~Diagnostics()
    {
    }
    ;
};

typedef boost::shared_ptr<Diagnostics > Diagnostics_sptr;

/// Diagnostics_basic provides the minimal set of statistical
/// quantities to be calculated for a Bunch.
class Diagnostics_basic : public Diagnostics
{
private:
    bool have_writers;
    std::string filename;
    Bunch_sptr bunch_sptr;
    Diagnostics_write_helper write_helper;
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

public:
    /// Create a Diagnostics_basic object
    /// @param bunch the Bunch
    /// @param filename filename for output
    Diagnostics_basic(Bunch_sptr bunch, std::string const& filename);

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_basic class is serial.
    virtual bool
    is_serial() const;

    virtual void
    init_writers(H5::H5File & file);

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
    get_bunchmin() const;

    virtual const MArray1d
    get_bunchmax() const;


    virtual void
    write();

    virtual
    ~Diagnostics_basic();
};

typedef boost::shared_ptr<Diagnostics_basic > Diagnostics_basic_sptr;

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_full2 : public Diagnostics
{
private:
    bool have_writers;
    std::string filename;
    Bunch_sptr bunch_sptr;
    Diagnostics_write_helper write_helper;
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
    Diagnostics_full2(Bunch_sptr bunch, std::string const& filename);

    virtual void
    init_writers(H5::H5File & file);

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

    virtual
    ~Diagnostics_full2();
};

typedef boost::shared_ptr<Diagnostics_full2 > Diagnostics_full2_sptr;

/// Diagnostics_particles dumps the state of particles in a bunch
class Diagnostics_particles : public Diagnostics
{
private:
    bool have_writers;
    H5::H5File file;
    int max_particles;
    Bunch_sptr bunch_sptr;
    std::string filename;
    Diagnostics_write_helper write_helper;
    void
    receive_other_local_particles(std::vector<int > const& local_nums, H5::H5File & file);
    void
    send_local_particles();
public:
    /// Create a Diagnostics_particles object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param max_particles the maximum number of particles to save (0 for all)
    Diagnostics_particles(Bunch_sptr bunch_sptr, std::string const& filename,
            int max_particles = 0);

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_particles class is not serial.
    virtual bool
    is_serial() const;

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();

    virtual
    ~Diagnostics_particles();
};

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
    Bunch_sptr bunch_sptr;
    std::string filename;
    Diagnostics_write_helper write_helper;
    double s;
    Hdf5_serial_writer<double > * writer_s;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double trajectory_length;
    Hdf5_serial_writer<double > * writer_trajectory_length;
    MArray1d coords;
    Hdf5_serial_writer<MArray1d_ref > * writer_coords;
    virtual void
    init_writers(H5::H5File & file);

public:
    /// Create an empty Diagnostics_track object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param particle_id the particle ID to track
    Diagnostics_track(Bunch_sptr bunch_sptr, std::string const& filename,
            int particle_id);

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_full2 class is serial.
    virtual bool
    is_serial() const;

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();

    virtual
    ~Diagnostics_track();
};

typedef boost::shared_ptr<Diagnostics_track > Diagnostics_track_sptr;

#endif /* DIAGNOSTICS_H_ */
