#ifndef DIAGNOSTICS_BASIC_H_
#define DIAGNOSTICS_BASIC_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_basic provides the minimal set of statistical
/// quantities to be calculated for a Bunch.
class Diagnostics_basic : public Diagnostics {
  public:
    static const char name[];

  private:
    bool have_writers;
    double s_n;
    Hdf5_serial_writer<double>* writer_s_n;
    int repetition;
    Hdf5_serial_writer<int>* writer_repetition;
    double s;
    Hdf5_serial_writer<double>* writer_s;
    int num_particles;
    Hdf5_serial_writer<int>* writer_num_particles;
    double real_num_particles;
    Hdf5_serial_writer<double>* writer_real_num_particles;
    double pz;
    Hdf5_serial_writer<double>* writer_pz;
    MArray1d mean;
    Hdf5_serial_writer<MArray1d_ref>* writer_mean;
    MArray1d std;
    Hdf5_serial_writer<MArray1d_ref>* writer_std;
    MArray1d min;
    Hdf5_serial_writer<MArray1d_ref>* writer_min;
    MArray1d max;
    Hdf5_serial_writer<MArray1d_ref>* writer_max;

  public:
    /// Create a Diagnostics_basic object
    /// @param filename filename for output
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_basic(std::string const& filename,
                      std::string const& local_dir = "");

    // Default constructor for serialization use only
    Diagnostics_basic();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_basic class is serial.
    virtual bool is_serial() const;

    virtual void init_writers(Hdf5_file_sptr file_sptr);

    /// Update the diagnostics
    virtual void update();

    /// Get the distance from the origin along the reference trajectory in
    /// meters.
    virtual double get_s_n() const;

    /// Get the number of complete repetitions.
    virtual int get_repetition() const;

    /// Get the total distance along the reference trajectory in meters.
    virtual double get_s() const;

    /// Get the total number of macroparticles in the bunch
    virtual int get_num_particles() const;

    /// Get the total number of real particles represented by the bunch
    virtual double get_real_num_particles() const;

    /// Get the current momentum of the bunch
    virtual double get_pz() const;

    /// Get a six-dimensional vector of the means of each phase-space
    /// coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref get_mean() const;

    /// Get a six-dimensional vector of the standard deviations of each
    /// phase-space coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref get_std() const;

    virtual const MArray1d get_min() const;

    virtual const MArray1d get_max() const;

    virtual void write();

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);

    virtual ~Diagnostics_basic();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_basic)
typedef boost::shared_ptr<Diagnostics_basic>
    Diagnostics_basic_sptr; // syndoc:include

#endif /* DIAGNOSTICS_BASIC_H_ */
