#ifndef DIAGNOSTICS_NF_H_
#define DIAGNOSTICS_NF_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/simulation/fast_normal_form.h"

/// Diagnostics_full2 provides the full set of statistical
/// quantities to be calculated for a Bunch up to the second moments.
class Diagnostics_normal_form : public Diagnostics
{
public:
    static const char name[];
private:
    Fast_normal_form_sptr fnf_sptr;
    bool have_writers;
    double s_n;
    Hdf5_serial_writer<double > * writer_s_n;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double s;
    Hdf5_serial_writer<double > * writer_s;
    int num_particles;
    Hdf5_serial_writer<int > * writer_num_particles;
    MArray1d aa2;
    Hdf5_serial_writer<MArray1d_ref > * writer_aa2;
    virtual void
    update_aa2();


public:
   
    Diagnostics_normal_form(Fast_normal_form_sptr fast_nf_sptr, std::string const& filename, std::string const& local_dir = "");

    // Default constructor for serialization use only
    Diagnostics_normal_form();

     virtual void
     init_writers(Hdf5_file_sptr file_sptr);

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_normal_form class is serial.
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
    ~Diagnostics_normal_form();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_normal_form)
typedef boost::shared_ptr<Diagnostics_normal_form > Diagnostics_normal_form_sptr; // syndoc:include

#endif /* DIAGNOSTICS_NF_H_ */
