#ifndef SPACE_CHARGE_RECTANGULAR_H_
#define SPACE_CHARGE_RECTANGULAR_H_

#include "synergia/simulation/operator.h"
#include "synergia/simulation/collective_operator_options.h"
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/utils/distributed_fft3d_rect.h"

//#include "synergia/collective/diagnostics_space_charge.h"
//#include  "synergia/collective/fftw_rectangular_helper.h"

// #include "synergia/utils/distributed_fft3d.h"


struct Space_charge_rectangular_options : public CO_options
{
    std::array<int, 3> shape;
    std::array<double, 3> pipe_size;
    int comm_group_size;

    Space_charge_rectangular_options(
            std::array<int, 3> const& shape = {32, 32, 64},
            std::array<double, 3> const& pipe_size = {0.1, 0.1, 1.0})
        : shape(shape), pipe_size(pipe_size), comm_group_size(1)
    { }

    CO_options * clone() const override
    { return new Space_charge_rectangular_options(*this); }

    Collective_operator * create_operator() const override;

    template<class Archive>
    void serialize(Archive & ar)
    {
        ar(cereal::base_class<CO_options>(this));
        ar(shape);
        ar(pipe_size);
        ar(comm_group_size);
    }
};

CEREAL_REGISTER_TYPE(Space_charge_rectangular_options)


class Space_charge_rectangular : public Collective_operator
{

private: 

    const Space_charge_rectangular_options options;

    Rectangular_grid_domain domain;

    Distributed_fft3d_rect fft;
    Commxx comm;

    karray1d_dev rho;
    karray1d_dev phi;
    karray1d_dev phihat;

#if 0
    karray1d_dev rho2hat;
    karray1d_dev phi2hat;
    karray1d_dev g2hat;
#endif

    karray1d_hst h_rho;
    karray1d_hst h_phi;

    karray1d_dev enx;
    karray1d_dev eny;
    karray1d_dev enz;

#if 0
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
    /// Example: consider  a bunch comunicator of size 128 and a spc comm optimal size =32 (good for the tev fermilab  cluster)
    /// when equally_spread=false, the spc comunicator will be on ranks [0,31]. ranks [32,127] do not have a spc comunicator
    /// -------------------------  the fourier transforms are done on the [0,31] ranks, and the rest of ranks wait doing nothing
    /// when equally_spread=true, there will be four spc comunicators, on ranks [0,31],[32,63],[64,95] and [96,127]
    ///------------------------- the fourier transforms are repeated independently on all spc comunicators     
    bool equally_spread;
    Diagnostics_space_charge_rectangulars diagnostics_list;
    bool have_diagnostics;
    
    
    void
    fill_guards_pplanes(Distributed_rectangular_grid & phi, int lower, int upper, int lengthx,
                          MArray2d & g_lower, MArray2d &g_upper);
    void
    construct_fftw_helper(Commxx_sptr comm_sptr);

#endif

private:

    void apply_impl( 
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger);

    void apply_bunch( 
            Bunch & bunch, 
            double time_step, 
            Logger & logger);

    void setup_communication(
            Commxx const & bunch_comm);

    void construct_workspaces(
            std::array<int, 3> const& s);

    void set_domain(
            Bunch const & bunch);

    void get_local_charge_density(
            Bunch const& bunch);

    void get_global_charge_density(
            Bunch const & bunch );

    void get_local_phi(double gamma);
    void get_global_phi();
    void extract_force();
    double get_normalization_force();

    void apply_kick(
            Bunch & bunch,
            double fn_norm,
            double time_step );



public:

    Space_charge_rectangular(
            Space_charge_rectangular_options const& ops);

#if 0
    Space_charge_rectangular(
            Commxx_sptr comm_f_sptr, 
            std::vector<double> const& pipe_size, 
            std::vector<int> const & grid_shape, 
            bool equally_spread);

    Space_charge_rectangular(
            std::vector<double> const& pipe_size, 
            std::vector<int > const & grid_shape);

    
    void
    set_fftw_helper(Commxx_sptr comm_sptr, bool equally_spread);
 
    bool 
    get_have_fftw_helper() const;

    
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
    
#endif
};

inline Collective_operator * 
Space_charge_rectangular_options::create_operator() const
{ return new Space_charge_rectangular(*this); }


#endif /* SPACE_CHARGE_RECTANGULAR_H_ */
