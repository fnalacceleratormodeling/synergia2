#ifndef DIAGNOSTICS_PHASE_SPACE_DENSITY_H_
#define DIAGNOSTICS_PHASE_SPACE_DENSITY_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_phase_space_density provides density and dipole density in
/// the longitudinal phase space by depositing the particles on a grid

class Diagnostics_phase_space_density : public Diagnostics
{
public:
    static const char name[];
private:
     int grid_z;
     int grid_zp;
     double z_nsigma;
     double zp_nsigma;
     double grid_zrms;
     double grid_zprms;
     double z_range;
     double zp_range;
     double z_cell;
     double zp_cell;     
     bool have_grid;
     bool have_writers;
     double s_n;
     Hdf5_serial_writer<double > * writer_s_n;
     int repetition;
     Hdf5_serial_writer<int > * writer_repetition;
     double s;     
     Hdf5_serial_writer<double > * writer_s;    
     MArray2d density;
     Hdf5_serial_writer<MArray2d_ref > * writer_density;
     MArray2d xdensity;
     Hdf5_serial_writer<MArray2d_ref > * writer_xdensity;
    // MArray2d ydensity;
     //Hdf5_serial_writer<MArray2d_ref > * writer_ydensity;



    virtual void
    init_writers(Hdf5_file_sptr file_sptr);
    void
    deposit_densities();

public:
    /// Create a Diagnostics_phase_space_density object
    /// @param bunch the Bunch
    /// @param filename filename for output
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_phase_space_density(std::string const& filename, int grid_z, int grid_zp, double  z_nsigma=4.0,
    double zp_nsigma=4.0, std::string const& local_dir = "");

    // Default constructor for serialization use only
    Diagnostics_phase_space_density();
    
    
    
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
    ~Diagnostics_phase_space_density();
};
 BOOST_CLASS_EXPORT_KEY(Diagnostics_phase_space_density)
typedef boost::shared_ptr<Diagnostics_phase_space_density > Diagnostics_phase_space_density_sptr; // syndoc:include

#endif /* DIAGNOSTICS_FULL2_H_ */
