#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include "components/bunch/bunch.h"
#include "utils/hdf5_writer.h"

/// Diagnostics provides the minimal set of statistical
/// quantities to be calculated for a Bunch.
class Diagnostics
{
private:
    bool have_writers;
protected:
    double s;
    Hdf5_writer<double > * writer_s;
    int repetition;
    Hdf5_writer<int > * writer_repetition;
    double trajectory_length;
    Hdf5_writer<double > * writer_trajectory_length;
    MArray1d mean;
    Hdf5_writer<MArray1d_ref > * writer_mean;
    MArray1d std;
    Hdf5_writer<MArray1d_ref > * writer_std;
    virtual void
    update_mean(Bunch const& bunch);
    virtual void
    update_std(Bunch const& bunch);
public:
    /// Create an empty Diagnostics object
    Diagnostics();

    /// Create a Diagnostics object
    /// @param bunch the Bunch
    Diagnostics(Bunch const& bunch);

    /// Update the diagnostics
    /// @param bunch the Bunch
    virtual void
    update(Bunch const& bunch);

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

    /// Get a six-dimensional vector of the means of each phase-space
    /// coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref
    get_mean() const;

    /// Get a six-dimensional vector of the standard deviations of each
    /// phase-space coordinate. The units are in Synergia units.
    virtual Const_MArray1d_ref
    get_std() const;

    virtual void
    init_writers(hid_t & hdf5_file);

    virtual void
    write_hdf5();

    virtual
    ~Diagnostics();
};

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_full2 : public Diagnostics
{
protected:
    MArray2d mom2;
    Hdf5_writer<MArray2d_ref > * writer_mom2;
    MArray2d corr;
    Hdf5_writer<MArray2d_ref > * writer_corr;
    double emitx, emity, emitz, emitxy, emitxyz;
    Hdf5_writer<double > *writer_emitx, *writer_emity, *writer_emitz,
            *writer_emitxy, *writer_emitxyz;
    virtual void
    update_full2(Bunch const& bunch);
    virtual void
    update_emittances();
public:
    /// Create an empty Diagnostics_full2 object
    Diagnostics_full2();

    /// Create a Diagnostics object
    /// @param bunch the Bunch
    Diagnostics_full2(Bunch const& bunch);

    /// Update the diagnostics
    /// @param bunch the Bunch
    virtual void
    update(Bunch const& bunch);

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
    init_writers(hid_t & hdf5_file);

    virtual void
    write_hdf5();

    virtual
    ~Diagnostics_full2();
};

#endif /* DIAGNOSTICS_H_ */
