#ifndef SPACE_CHARGE_RECTANGULAR_H_
#define SPACE_CHARGE_RECTANGULAR_H_
// #include "synergia/simulation/operator.h"
#include "synergia/collective/diagnostics_space_charge.h"
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
    ///pipe size, x,y,z meters, lab frame
    std::vector<double > pipe_size; 
    std::vector<int > grid_shape;
    Rectangular_grid_domain_sptr domain_sptr;
    /// communicator for parallel Fourier transforms 
    Commxx_sptr comm_f_sptr;
    Fftw_rectangular_helper_sptr fftw_helper_sptr;
    bool have_fftw_helper;
    bool have_domain;
    /// the spc communicator is constructed  on a subset of bunch communicator ranks, with an optimal size
    /// usually proportional to the number of processors of a node
    /// Example: consider  a bunch communicator of size 128 and a spc comm optimal size =32 (good for the tev fermilab  cluster)
    /// when equally_spread=false, the spc communicator will be on ranks [0,31]. ranks [32,127] do not have a spc communicator
    /// -------------------------  the fourier transforms are done on the [0,31] ranks, and the rest of ranks wait doing nothing
    /// when equally_spread=true, there will be four spc communicators, on ranks [0,31],[32,63],[64,95] and [96,127]
    ///------------------------- the fourier transforms are repeated independently on all spc communicators
    bool equally_spread;
    Diagnostics_space_charge_rectangulars diagnostics_list;
    bool have_diagnostics;
    
    
    void
    fill_guards_pplanes(Distributed_rectangular_grid & phi, int lower, int upper, int lengthx,
                          MArray2d & g_lower, MArray2d &g_upper);
    void
    construct_fftw_helper(Commxx_sptr comm_sptr);
    

public:
    Space_charge_rectangular(Commxx_sptr comm_f_sptr, std::vector<double > const & pipe_size, 
			       std::vector<int > const & grid_shape, bool equally_spread);
    Space_charge_rectangular(std::vector<double > const & pipe_size, std::vector<int > const & grid_shape);
    Space_charge_rectangular();
    
    virtual Space_charge_rectangular *
    clone();

    void
    set_fftw_helper(Commxx_sptr comm_sptr, bool equally_spread);
 
    bool 
    get_have_fftw_helper() const;

    void 
    set_domain(Bunch const & bunch);
    
    bool
    get_have_domain() const;


    Commxx_sptr 
    get_comm_sptr() const;
    
    
    std::vector<double >
    get_pipe_size() const;

    std::vector<int >
    get_grid_shape() const;
    
    
    void
    set_diagnostics_list(Diagnostics_space_charge_rectangulars diagnosticss);
   
    void
    add_diagnostics(Diagnostics_space_charge_rectangular_sptr diagnostics_sptr);

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
    get_phi_local(Rectangular_grid & rho, double const& gamma);

    Rectangular_grid_sptr
    get_En( Distributed_rectangular_grid & phi_local, int component);

    std::vector<Rectangular_grid_sptr>
    get_Efield(Rectangular_grid & rho,Bunch const& bunch, int max_component, double const & gamma);
      
    void
    do_diagnostics(Rectangular_grid const& En, int component, double time_step, Step & step, Bunch & bunch);
    
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

typedef boost::shared_ptr<Space_charge_rectangular> Space_charge_rectangular_sptr; // syndoc:include

#endif /* SPACE_CHARGE_RECTANGULAR_H_ */
