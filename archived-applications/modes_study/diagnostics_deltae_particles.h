#ifndef DIAGNOSTICS_DELTAE_PARTICLES_H_
#define DIAGNOSTICS_DELTAE_PARTICLES_H_

#include "synergia/bunch/diagnostics.h"

class Diagnostics_deltae_particles : public Diagnostics
{
  public:
    static const char name[];
  private:
    int turn;
    int begin_turn;
    int end_turn;
    int average_bin;
    double beta_x;
    double alpha_x;
    double beta_y;
    double alpha_y;
    double beta_z;
    double alpha_z;
    bool have_writers;
    bool have_initial_deltae;
    MArray2d *local_deltae;
  

    void    
    send_local_particles();

    void
    receive_other_local_particles(std::vector<int > const& local_nums,
            Hdf5_file_sptr file_sptr);
   
    std::string
    get_deltae_serialization_path() const;
    
  public:
     Diagnostics_deltae_particles(std::string const& filename, int begin_turn, int end_turn,  int average_bin, 
                               double beta_x, double alpha_x, double beta_y, double alpha_y, 
                                 double beta_z, double alpha_z, std::string const& local_dir = "");
     
     Diagnostics_deltae_particles();
     
    virtual bool
    is_serial() const;

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();
    
   
    
   
   template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const;

   template<class Archive>
        void
        load(Archive & ar, const unsigned int version);
   BOOST_SERIALIZATION_SPLIT_MEMBER() 

   virtual
    ~Diagnostics_deltae_particles();
   
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_deltae_particles)
typedef boost::shared_ptr<Diagnostics_deltae_particles > Diagnostics_deltae_particles_sptr; // syndoc:include

#endif /* DIAGNOSTICS_DELTAE_PARTICLES_H_ */