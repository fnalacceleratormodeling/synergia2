#ifndef SPACE_CHARGE_2D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_2D_OPEN_HOCKNEY_H_

#include "synergia/simulation/operator.h"
#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/utils/distributed_fft2d.h"


struct Space_charge_2d_open_hockney_options : public CO_options
{
    std::array<int, 3> shape;
    std::array<int, 3> doubled_shape;

    bool periodic_z;
    double z_period;
    bool grid_entire_period;
    double n_sigma;

    int comm_group_size;

    Space_charge_2d_open_hockney_options(int gridx, int gridy, int gridz)
        : shape{gridx, gridy, gridz}
        , doubled_shape{gridx*2, gridy*2, gridz}
        , periodic_z(false)
        , z_period(0.0)
        , grid_entire_period(false)
        , n_sigma(8.0)
        , comm_group_size(4)
    { }

    Collective_operator * create_operator() const override;
};


class Space_charge_2d_open_hockney : public Collective_operator
{

public:

#if 0
    enum Green_fn_type
    {
        pointlike = 1, bruteforce = 2
    };
#endif

private:

    Space_charge_2d_open_hockney_options options;

    Rectangular_grid_domain domain;
    Rectangular_grid_domain doubled_domain;

    karray2d_dev particle_bin;

    Distributed_fft2d fft;
    Commxx comm;

#if 0
    std::vector<int > grid_shape, doubled_grid_shape;
    Rectangular_grid_domain_sptr domain_sptr, doubled_domain_sptr;
    bool periodic_z;
    double z_period;
    bool grid_entire_period;
    Green_fn_type green_fn_type;
    Distributed_fft2d_sptr distributed_fft2d_sptr;

    Commxx_divider_sptr commxx_divider_sptr;
    Commxx_sptr comm2_sptr, comm1_sptr;

    std::vector<int > lowers1, lengths1;
    int real_lower, real_upper, real_length;
    std::vector<int > real_lengths;
    std::vector<int > real_lengths_1d;

    int doubled_lower, doubled_upper;
    int real_doubled_lower, real_doubled_upper;

    double n_sigma;
    double bunch_particle_charge, bunch_total_num;
    double beta, gamma;
    bool use_cell_coords;
    bool need_state_conversion;
    bool domain_fixed;
    bool have_domains;

    void setup_communication(Commxx_sptr const& bunch_comm_sptr);
    void setup_derived_communication();
    void setup_default_options();
    void set_doubled_domain();
#endif

    void apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger) override;

    void apply_bunch(
            Bunch & bunch, 
            double time_step, 
            Logger & logger);

    void setup_communication(Commxx const & bunch_comm);
    void update_domain(Bunch const & bunch);

    karray1d_dev get_local_charge_density(Bunch const& bunch);
    karray1d_dev get_green_fn2_pointlike();
    karray1d_dev get_local_force2(karray1d_dev & rho2, karray1d_dev & g2);

    void get_global_charge_density(karray1d_dev & rho2, Bunch const & bunch);
    void get_global_force2(karray1d_dev & phi2);

    void apply_kick(
            Bunch & bunch,
            karray1d_dev const & rho2,
            karray1d_dev const & fn2,
            double fn_norm,
            double time_step );

    double get_normalization_force(Bunch const & bunch);

public:

    Space_charge_2d_open_hockney(Space_charge_2d_open_hockney_options const & ops);

#if 0
    bool get_need_state_conversion();
    double get_n_sigma() const;

    void set_fixed_domain(Rectangular_grid_domain_sptr domain_sptr);
    void update_domain(Bunch const& bunch);

    Rectangular_grid_domain const& get_domain() const
    { return *domain_sptr; }

    /// Returns local charge density on doubled grid in [C/m^3]
    Rectangular_grid_sptr
    get_local_charge_density(Bunch const& bunch);

    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2(
            Rectangular_grid const& local_charge_density,
            Commxx_sptr comm_sptr);

    /// Returns Green function on the doubled grid in [1/m^3]
    Distributed_rectangular_grid_sptr
    get_green_fn2_pointlike();

    Distributed_rectangular_grid_sptr
    get_green_fn2_brute_force();

    Distributed_rectangular_grid_sptr
    get_local_force2(Distributed_rectangular_grid & charge_density2,
            Distributed_rectangular_grid & green_fn2);

    Rectangular_grid_sptr
    get_global_electric_force2(
            Distributed_rectangular_grid const& dist_force);

    void
    apply_kick(Bunch & bunch, Distributed_rectangular_grid const& rho2_1d,
            Rectangular_grid const& Fn, double delta_tau);

    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
#endif
};

inline Collective_operator * 
Space_charge_2d_open_hockney_options::create_operator() const
{ return new Space_charge_2d_open_hockney(*this); }

#endif /* SPACE_CHARGE_2D_OPEN_HOCKNEY_H_ */
