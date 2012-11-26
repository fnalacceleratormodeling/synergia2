#ifndef SPACE_CHARGE_3D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_3D_OPEN_HOCKNEY_H_
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/distributed_rectangular_grid.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/distributed_fft3d.h"

/// Note: internal grid is stored in [z][y][x] order, but
/// grid shape expects [x][y][z] order.
class Space_charge_3d_open_hockney : public Collective_operator
{
public:
    enum Green_fn_type
    {
        pointlike = 1, linear = 2
    };
    enum Charge_density_comm
    {
        reduce_scatter = 1, charge_allreduce = 2
    };
    enum E_field_comm
    {
        gatherv_bcast = 1, allgatherv = 2, e_field_allreduce = 3
    };
private:
    std::vector<int > grid_shape, doubled_grid_shape, padded_grid_shape;
    Rectangular_grid_domain_sptr domain_sptr, doubled_domain_sptr;
    bool periodic_z;
    double z_period;
    bool grid_entire_period;
    bool longitudinal_kicks;
    Green_fn_type green_fn_type;
    Charge_density_comm charge_density_comm;
    E_field_comm e_field_comm;
    Distributed_fft3d_sptr distributed_fft3d_sptr;
    Commxx_sptr comm2_sptr, comm1_sptr;
    std::vector<int > lowers1, lengths1;
    int real_lower, real_upper, real_length;
    std::vector<int > real_lengths;
    int doubled_lower, doubled_upper;
    int real_doubled_lower, real_doubled_upper;
    double n_sigma;
    bool domain_fixed;
    bool have_domains;
    void
    setup_nondoubled_communication();
    void
    setup_default_options();
    void
    set_doubled_domain();
public:
    Space_charge_3d_open_hockney(Commxx_sptr comm_sptr,
            std::vector<int > const & grid_shape,
            bool longitudinal_kicks = true, bool periodic_z = false,
            double z_period = 0.0, bool grid_entire_period = false,
            double n_sigma = 8.0);
    /// Note: Use Space_charge_3d_open_hockney::get_internal_grid_shape for
    /// Distributed_fft3d.
    Space_charge_3d_open_hockney(Distributed_fft3d_sptr distributed_fft3d_sptr,
            bool longitudinal_kicks = true, bool periodic_z = false,
            double z_period = 0.0, bool grid_entire_period = false,
            double n_sigma = 8.0);
    Space_charge_3d_open_hockney();
    double
    get_n_sigma() const;
    void
    set_green_fn_type(Green_fn_type green_fn_type);
    Green_fn_type
    get_green_fn_type() const;
    void
    set_charge_density_comm(Charge_density_comm charge_density_comm);
    Charge_density_comm
    get_charge_density_comm() const;
    void
    set_e_field_comm(E_field_comm e_field_comm);
    E_field_comm
    get_e_field_comm() const;
    void
    auto_tune_comm(bool verbose = false);
    void
    set_fixed_domain(Rectangular_grid_domain_sptr domain_sptr);
    void
    update_domain(Bunch const& bunch);
    Rectangular_grid_domain const&
    get_domain() const
    {
        return *domain_sptr;
    }
    Rectangular_grid_domain_sptr
    get_domain_sptr() const;
    Rectangular_grid_domain_sptr
    get_doubled_domain_sptr() const;
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2_reduce_scatter(
            Rectangular_grid const& local_charge_density, Commxx_sptr comm_sptr);
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2_allreduce(
            Rectangular_grid const& local_charge_density, Commxx_sptr comm_sptr);
    /// Returns local charge density on original grid in [C/m^3]
    Rectangular_grid_sptr
    get_local_charge_density(Bunch const& bunch);
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2(Rectangular_grid const& local_charge_density,
            Commxx_sptr comm_sptr);
    /// Returns Green function on the doubled grid in [1/m^3]
    Distributed_rectangular_grid_sptr
    get_green_fn2_pointlike();
    Distributed_rectangular_grid_sptr
    get_green_fn2_linear();
    Distributed_rectangular_grid_sptr
    get_scalar_field2(Distributed_rectangular_grid & charge_density22,
            Distributed_rectangular_grid & green_fn2);
    Distributed_rectangular_grid_sptr
    extract_scalar_field(Distributed_rectangular_grid const& scalar_field2);
    /// Returns component of electric field [V/m]
    /// @param scalar_field the scalar field [V]
    /// @param component which component (0=x, 1=y, 2=z)
    Distributed_rectangular_grid_sptr
    get_electric_field_component(
            Distributed_rectangular_grid const& scalar_field, int component);
    Rectangular_grid_sptr
    get_global_electric_field_component_gatherv_bcast(
            Distributed_rectangular_grid const& dist_field);
    Rectangular_grid_sptr
    get_global_electric_field_component_allgatherv(
            Distributed_rectangular_grid const& dist_field);
    Rectangular_grid_sptr
    get_global_electric_field_component_allreduce(
            Distributed_rectangular_grid const& dist_field);
    Rectangular_grid_sptr
    get_global_electric_field_component(
            Distributed_rectangular_grid const& dist_field);
    void
    apply_kick(Bunch & bunch, Rectangular_grid const& En, double delta_tau,
            int component);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger);
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const;
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version);
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    virtual
    ~Space_charge_3d_open_hockney();
};
BOOST_CLASS_EXPORT_KEY(Space_charge_3d_open_hockney)

typedef boost::shared_ptr<Space_charge_3d_open_hockney > Space_charge_3d_open_hockney_sptr; // syndoc:include

#endif /* SPACE_CHARGE_3D_OPEN_HOCKNEY_H_ */
