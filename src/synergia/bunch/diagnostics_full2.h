#ifndef DIAGNOSTICS_FULL2_H_
#define DIAGNOSTICS_FULL2_H_

#include "synergia/bunch/diagnostics.h"
#include "boost/multi_array.hpp"

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_full2 : public Diagnostics
{
public:
    static const char name[];
private:
    bool have_writers = false;
    double s_n = 0.0;
    Hdf5_serial_writer<double > * writer_s_n = nullptr;
    int repetition = 0;
    Hdf5_serial_writer<int > * writer_repetition = nullptr;
    double s = 0.0;
    Hdf5_serial_writer<double > * writer_s = nullptr;
    int num_particles = 0;
    Hdf5_serial_writer<int > * writer_num_particles = nullptr;
    double real_num_particles = 0.0;
    Hdf5_serial_writer<double > * writer_real_num_particles = nullptr;
    double pz = 0.0;
    Hdf5_serial_writer<double > * writer_pz = nullptr;
    MArray1d mean {boost::extents[6]};
    Hdf5_serial_writer<MArray1d_ref > * writer_mean = nullptr;
    MArray1d std {boost::extents[6]};
    Hdf5_serial_writer<MArray1d_ref > * writer_std = nullptr;
    MArray1d min {boost::extents[3]};
    Hdf5_serial_writer<MArray1d_ref > * writer_min = nullptr;
    MArray1d max {boost::extents[3]};
    Hdf5_serial_writer<MArray1d_ref > * writer_max = nullptr;
    MArray2d mom2 {boost::extents[6][6]};
    Hdf5_serial_writer<MArray2d_ref > * writer_mom2 = nullptr;
    MArray2d corr {boost::extents[6][6]};
    Hdf5_serial_writer<MArray2d_ref > * writer_corr = nullptr;
    double emitx = 0.0;
    double emity = 0.0;
    double emitz = 0.0;
    double emitxy = 0.0;
    double emitxyz = 0.0;
    Hdf5_serial_writer<double >* writer_emitx = nullptr;
    Hdf5_serial_writer<double>* writer_emity = nullptr;
    Hdf5_serial_writer<double>* writer_emitz = nullptr;
    Hdf5_serial_writer<double>* writer_emitxy = nullptr;
    Hdf5_serial_writer<double>* writer_emitxyz = nullptr;
    
    virtual void update_full2();
    virtual void update_emittances();

public:
    /// Create a Diagnostics_full2 object
    /// @param bunch the Bunch
    /// @param filename filename for output
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_full2(std::string const& filename, std::string const& local_dir = "");

    // Default constructor for serialization use only
    Diagnostics_full2();

    Diagnostics_full2(Diagnostics_full2 const&) = delete;
    Diagnostics_full2& operator=(Diagnostics_full2 const&) = delete;

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
    get_s_n() const;

    /// Get the number of complete repetitions.
    virtual int
    get_repetition() const;

    /// Get the total distance along the reference trajectory in meters.
    virtual double
    get_s() const;

    /// Get the total number of macroparticles in the bunch
    virtual int
    get_num_particles() const;

    /// Get the total number of real particles represented by the bunch
    virtual double
    get_real_num_particles() const;

    /// Get the current momentum of the bunch
    virtual double
    get_pz() const;

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
typedef boost::shared_ptr<Diagnostics_full2 > Diagnostics_full2_sptr; // syndoc:include

#endif /* DIAGNOSTICS_FULL2_H_ */
