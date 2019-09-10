#ifndef DIAGNOSTICS_LOSS_H_
#define DIAGNOSTICS_LOSS_H_

#include "synergia/bunch/diagnostics.h"


class Diagnostics_loss : public Diagnostics
{ 
public:  

    constexpr static const char* diag_type = "diagnostics_loss";
    constexpr static const bool  diag_write_serial = true;

private:

    int bucket_index;
    int repetition;
    double s_ref_particle;
    double sn_ref_particle;

    karray2d_row coords;
    
private:

     void do_update(Bunch const& bunch) override {};
     void do_write (Bunch const& bunch) override;
    
public: 

    Diagnostics_loss( 
            std::string const& filename, 
            std::string const& local_dir = "" );
    
    void update(Bunch const& bunch, karray2d_row coords);
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#endif /* LOSS_DIAGNOSTICS_H_ */
 
