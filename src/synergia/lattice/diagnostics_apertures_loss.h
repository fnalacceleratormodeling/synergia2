#ifndef APLOSS_DIAGNOSTICS_H_
#define APLOSS_DIAGNOSTICS_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/multi_array_typedefs.h"
#include <boost/shared_ptr.hpp>
#include "synergia/utils/hdf5_file.h"


class Diagnostics_apertures_loss : public Diagnostics
{ 
    static const char name[];
private:
    // Commxx_sptr comm_sptr;     
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
    Diagnostics_apertures_loss(std::string const& filename, std::string const& local_dir = "");
    Diagnostics_apertures_loss();
    
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

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
   

     virtual
    ~Diagnostics_apertures_loss();
    
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_apertures_loss)
typedef boost::shared_ptr< Diagnostics_apertures_loss> Diagnostics_apertures_loss_sptr;
typedef std::list<Diagnostics_apertures_loss_sptr >  Diagnostics_apertures_losses;


 #endif /* APLOSS_DIAGNOSTICS_H_ */
 
