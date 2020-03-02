#ifndef DIAGNOSTICS_LOSS_H_
#define DIAGNOSTICS_LOSS_H_

#include "synergia/bunch/diagnostics.h"


class Diagnostics_loss : public Diagnostics
{
private:

    int bucket_index;
    int repetition;
    double s_ref_particle;
    double sn_ref_particle;

    karray2d_row coords;
 
public:

    Diagnostics_loss()
        : Diagnostics("diagnostics_loss", true)
    { }

private:

    void do_update(Bunch const& bunch) override;
    void do_reduce(Commxx const& comm, int root) override;
    void do_first_write(Hdf5_file& file) override;
    void do_write(Hdf5_file& file) override;

    friend class cereal::access;

    template<class Archive>
    void serialize(Archive & ar)
    {
        ar(cereal::base_class<Diagnostics>(this));
        ar(bucket_index);
    }
};

CEREAL_REGISTER_TYPE(Diagnostics_loss)

#if 0
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
    
    std::unique_ptr<Diagnostics> do_pilfer() override
    { return std::make_unique<Diagnostics_loss>(std::move(*this)); }

public: 

    Diagnostics_loss( 
            std::string const& filename = "", 
            std::string const& local_dir = "" );
    
    void update(Bunch const& bunch, karray2d_row coords);
    
    template<class Archive>
    void serialize(Archive & ar)
    {
        ar(cereal::base_class<Diagnostics>(this));
        ar(bucket_index);
    }
};
#endif

//CEREAL_REGISTER_TYPE(Diagnostics_loss)

#endif /* LOSS_DIAGNOSTICS_H_ */
 
