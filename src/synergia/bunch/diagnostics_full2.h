#ifndef DIAGNOSTICS_FULL2_H_
#define DIAGNOSTICS_FULL2_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_full2 : public Diagnostics
{
public:
    static const char name[];
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
    /// Create a Diagnostics_full2 object
    /// @param bunch the Bunch
    /// @param filename filename for output
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_full2(std::string const& filename, std::string const& local_dir = "");

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

    //begin egs screwing around
    bool get_have_writers();
    // end egs screwing around
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Diagnostics_full2();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_full2)
typedef boost::shared_ptr<Diagnostics_full2 > Diagnostics_full2_sptr;

#endif /* DIAGNOSTICS_FULL2_H_ */
