#ifndef LOSS_DIAGNOSTICS_H_
#define LOSS_DIAGNOSTICS_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/multi_array_typedefs.h"
#include <boost/shared_ptr.hpp>
#include "synergia/utils/hdf5_file.h"


class Diagnostics_loss : public Diagnostics
{ 
public:  
    static const char name[];
    static const char aperture_type[];//="aperture"
    static const char zcut_type[]; // ="zcut"
private:

    /// type choices are "zcut" or "ape"
     std::string type;  
     bool have_writers;
     std::vector<int> bucket_index;
     std::vector<int> repetition;
     std::vector<double> s_ref_particle;
     std::vector<double> sn_ref_particle;
     std::vector<MArray1d> coords;
     Hdf5_serial_writer<int > * writer_bucket_index;
     Hdf5_serial_writer<int > * writer_repetition;
     Hdf5_serial_writer<double > * writer_s;
     Hdf5_serial_writer<double > * writer_s_n;
     Hdf5_serial_writer<MArray1d_ref > * writer_coords;
    
 
    
public: 
    Diagnostics_loss(std::string const& filename, std::string const& type, std::string const& local_dir = "");
    Diagnostics_loss();
    
    bool
    is_serial() const;
    
    virtual void
    update( );
    
    virtual void
    update(int index, int rep, double s, double s_n,  MArray1d_ref coords );
  
    
    virtual void
    write();
    
    void
    init_writers(Hdf5_file_sptr file_sptr);

     std::string
     get_type() const;
    
    

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
   

     virtual
    ~Diagnostics_loss();
    
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_loss)
typedef boost::shared_ptr< Diagnostics_loss> Diagnostics_loss_sptr;
typedef std::list<Diagnostics_loss_sptr >  Diagnostics_losses;


 #endif /* LOSS_DIAGNOSTICS_H_ */
 
