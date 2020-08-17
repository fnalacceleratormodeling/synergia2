#ifndef SPACE_CHARGE_2D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_2D_OPEN_HOCKNEY_H_

#include "synergia/simulation/operator.h"
#include "synergia/simulation/collective_operator_options.h"

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

    Space_charge_2d_open_hockney_options(
            int gridx = 32, int gridy = 32, int gridz = 32)
        : shape{gridx, gridy, gridz}
        , doubled_shape{gridx*2, gridy*2, gridz}
        , periodic_z(false)
        , z_period(0.0)
        , grid_entire_period(false)
        , n_sigma(8.0)
        , comm_group_size(4)
    { }

    CO_options * clone() const override
    { return new Space_charge_2d_open_hockney_options(*this); }

    Collective_operator * create_operator() const override;

    template<class Archive>
    void serialize(Archive & ar)
    { 
        ar(cereal::base_class<CO_options>(this));
        ar(shape);
        ar(doubled_shape);
        ar(periodic_z);
        ar(z_period);
        ar(grid_entire_period);
        ar(n_sigma);
        ar(comm_group_size);
    }
};

CEREAL_REGISTER_TYPE(Space_charge_2d_open_hockney_options)


class Space_charge_2d_open_hockney : public Collective_operator
{

private:

    const Space_charge_2d_open_hockney_options options;

    Rectangular_grid_domain domain;
    Rectangular_grid_domain doubled_domain;

    karray2d_dev particle_bin;

    Distributed_fft2d fft;
    Commxx comm;

    karray1d_dev rho2;
    karray1d_dev phi2;
    karray1d_dev g2;

    karray1d_hst h_rho2;
    karray1d_hst h_phi2;

    void construct_workspaces(std::array<int, 3> const& s);

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

    void get_local_charge_density(Bunch const& bunch);
    void get_green_fn2_pointlike();
    void get_local_force2();

    void get_global_charge_density(Bunch const & bunch);
    void get_global_force2();

    void apply_kick(
            Bunch & bunch,
            double fn_norm,
            double time_step );

    double get_normalization_force(Bunch const & bunch);

public:

    Space_charge_2d_open_hockney(Space_charge_2d_open_hockney_options const & ops);

};

inline Collective_operator * 
Space_charge_2d_open_hockney_options::create_operator() const
{ return new Space_charge_2d_open_hockney(*this); }

#endif /* SPACE_CHARGE_2D_OPEN_HOCKNEY_H_ */
