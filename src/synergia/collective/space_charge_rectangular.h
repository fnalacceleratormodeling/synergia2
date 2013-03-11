#ifndef SPACE_CHARGE_RECTANGULAR_H_
#define SPACE_CHARGE_RECTANGULAR_H_
// #include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/distributed_rectangular_grid.h"
#include "synergia/utils/commxx.h"
#include  "synergia/collective/fftw_rectangular_helper.h"

// #include "synergia/utils/distributed_fft3d.h"


class Space_charge_rectangular : public Collective_operator
{
private:
    std::vector<double > pipe_size; //pipe size, x,y,x meters
    std::vector<int > grid_shape;
    Rectangular_grid_domain_sptr domain_sptr;
    Commxx_sptr comm_f_sptr;
    Fftw_rectangular_helper_sptr fftw_helper_sptr;
    bool have_fftw_helper;
    void
    fill_guards_pplanes(Distributed_rectangular_grid & phi, int lower, int upper, int lengthx,
                          MArray2d & g_lower, MArray2d &g_upper);
    void
    construct_fftw_helper(Commxx_sptr comm_sptr);

public:

    Space_charge_rectangular(Commxx_sptr comm_f_sptr, std::vector<double > const & pipe_size, std::vector<int > const & grid_shape);
    Space_charge_rectangular(std::vector<double > const & pipe_size, std::vector<int > const & grid_shape);
    Space_charge_rectangular();
    
    virtual Space_charge_rectangular *
    clone();

    void
    set_fftw_helper(Commxx_sptr comm_sptr);
    
    Commxx_sptr 
    get_comm_sptr() const;
    
    
    std::vector<double >
    get_pipe_size() const;

    std::vector<int >
    get_grid_shape() const;

   Rectangular_grid_domain const&
   get_domain() const
   {
       return *domain_sptr;
   }

   Rectangular_grid_domain_sptr
   get_domain_sptr() const;

   Fftw_rectangular_helper_sptr
   get_fftw_helper_sptr() const;


   Rectangular_grid_sptr
   get_charge_density(Bunch const& bunch);


   Distributed_rectangular_grid_sptr
   get_phi_local(Rectangular_grid & rho);

   Rectangular_grid_sptr
   get_En( Distributed_rectangular_grid & phi_local, int component);

   std::vector<Rectangular_grid_sptr>
   get_Efield(Rectangular_grid & rho,Bunch const& bunch, int max_component);
    
   void
   apply_kick(Bunch & bunch, Rectangular_grid const& En, double time_step, int component);

    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const;
   template<class Archive>
        void
        load(Archive & ar, const unsigned int version);            
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
    virtual
    ~Space_charge_rectangular();
};
BOOST_CLASS_EXPORT_KEY(Space_charge_rectangular);
typedef boost::shared_ptr<Space_charge_rectangular> Space_charge_rectangular_sptr;

#endif /* SPACE_CHARGE_RECTANGULAR_H_ */
