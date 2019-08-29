#ifndef DIAGNOSTICS_LOSS_H_
#define DIAGNOSTICS_LOSS_H_

#include "synergia/bunch/diagnostics.h"


class Diagnostics_loss : public Diagnostics
{ 
public:  

    constexpr static const char* diag_type = "diagnostics_loss";
    constexpr static const bool  diag_write_serial = true;

private:

    std::vector<int> bucket_index;
    std::vector<int> repetition;
    std::vector<double> s_ref_particle;
    std::vector<double> sn_ref_particle;
    std::vector<std::array<double, 7>> coords;
    
private:

     void do_update() override {};
     void do_write()  override;
    
public: 

    Diagnostics_loss( 
            std::string const& filename, 
            std::string const& local_dir = "" );
    
    void update(
            int index, int rep, double s, double s_n,  
            std::array<double, 7> const& coords );
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#endif /* LOSS_DIAGNOSTICS_H_ */
 
