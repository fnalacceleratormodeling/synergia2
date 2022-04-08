#ifndef SPACE_CHARGE_3D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_3D_OPEN_HOCKNEY_H_

#include "synergia/simulation/collective_operator_options.h"
#include "synergia/simulation/operator.h"

#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/rectangular_grid_domain.h"

#include "synergia/utils/distributed_fft3d.h"

enum class green_fn_t {
  pointlike,
  linear,
};

class Space_charge_3d_open_hockney;

struct Space_charge_3d_open_hockney_options
  : public CO_base_options<Space_charge_3d_open_hockney_options,
                           Space_charge_3d_open_hockney> {
  std::array<int, 3> shape;
  std::array<int, 3> doubled_shape;

  green_fn_t green_fn;
  bool periodic_z;
  double z_period;
  bool grid_entire_period;
  double n_sigma;
  double kick_scale;
  bool domain_fixed;

  int comm_group_size;

  Space_charge_3d_open_hockney_options(int gridx = 32,
                                       int gridy = 32,
                                       int gridz = 64)
    : shape{gridx, gridy, gridz}
    , doubled_shape{gridx * 2, gridy * 2, gridz * 2}
    , green_fn(green_fn_t::linear)
    , periodic_z(false)
    , z_period(0.0)
    , grid_entire_period(false)
    , n_sigma(8.0)
    , kick_scale(1.0)
    , domain_fixed(false)
    , comm_group_size(4)
  {}

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(cereal::base_class<CO_base_options>(this));
    ar(shape);
    ar(doubled_shape);
    ar(green_fn);
    ar(periodic_z);
    ar(z_period);
    ar(grid_entire_period);
    ar(n_sigma);
    ar(comm_group_size);
  }
};

CEREAL_REGISTER_TYPE(Space_charge_3d_open_hockney_options)

/// Note: internal grid is stored in [z][y][x] order, but
/// grid shape expects [x][y][z] order.
class Space_charge_3d_open_hockney : public Collective_operator {
private:
  const Space_charge_3d_open_hockney_options options;

  std::string bunch_sim_id;

  Rectangular_grid_domain domain;
  Rectangular_grid_domain doubled_domain;
  bool use_fixed_domain;

  std::array<std::vector<Distributed_fft3d>, 2> ffts;

  karray1d_dev rho2;
  karray1d_dev phi2;
  karray1d_dev g2;

  karray1d_hst h_rho2;
  karray1d_hst h_phi2;

  karray1d_dev enx;
  karray1d_dev eny;
  karray1d_dev enz;

private:
  void apply_impl(Bunch_simulator& simulator, double time_step, Logger& logger);

  void apply_bunch(Bunch& bunch,
                   Distributed_fft3d& fft,
                   double time_step,
                   Logger& logger);

  void construct_workspaces(Bunch_simulator const& sim);

  void update_domain(Bunch const& bunch);

  void get_local_charge_density(Bunch const& bunch);

  void get_global_charge_density(Bunch const& bunch);

  void apply_kick(Bunch& bunch, double fn_norm, double time_step);

  void get_green_fn2_pointlike();
  void get_green_fn2_linear();

  void get_local_phi2(Distributed_fft3d& fft);

  void get_global_phi2(Distributed_fft3d const& fft);

  void get_force();

  double get_normalization_force(Distributed_fft3d const& fft);

public:
  Space_charge_3d_open_hockney(Space_charge_3d_open_hockney_options const& ops);

  void set_fixed_domain(std::array<double, 3> offset,
                        std::array<double, 3> size);

#if 0
private:

    std::vector<int > grid_shape, doubled_grid_shape, padded_grid_shape;
    Rectangular_grid_domain_sptr domain_sptr, doubled_domain_sptr;
    bool periodic_z;
    double z_period;
    bool grid_entire_period;
    bool longitudinal_kicks;
    Green_fn_type green_fn_type;
    Charge_density_comm charge_density_comm;
    E_field_comm e_field_comm; // vestigial from when we thought we would choose the most
                               // most performant collective primitive
    Distributed_fft3d_sptr distributed_fft3d_sptr;
    Commxx_divider_sptr commxx_divider_sptr;
    Commxx_sptr comm2_sptr, comm1_sptr;
    std::vector<int > lowers1, lengths1;
    int real_lower, real_upper, real_length;
    std::vector<int > real_lengths;
    int doubled_lower, doubled_upper;
    int real_doubled_lower, real_doubled_upper;
    double n_sigma;
    bool domain_fixed;
    bool have_domains;
    Diagnostics_space_charge_3d_hockneys diagnostics_list;
    bool have_diagnostics;
    double kick_scale;
    void
    constructor_common(std::vector<int > const& grid_shape);
    void
    set_doubled_domain();

public:

    Space_charge_3d_open_hockney(std::vector<int > const & grid_shape,
            bool longitudinal_kicks = true, bool periodic_z = false,
            double z_period = 0.0, bool grid_entire_period = false,
            double n_sigma = 8.0, double kick_scale=1.0);
    Space_charge_3d_open_hockney(Commxx_divider_sptr commxx_divider_sptr,
            std::vector<int > const & grid_shape,
            bool longitudinal_kicks = true, bool periodic_z = false,
            double z_period = 0.0, bool grid_entire_period = false,
            double n_sigma = 8.0, double kick_scale=1.0);
    /// Deprecated. The comm_sptr argument is ignored.
    Space_charge_3d_open_hockney(Commxx_sptr comm_sptr,
            std::vector<int > const & grid_shape,
            bool longitudinal_kicks = true, bool periodic_z = false,
            double z_period = 0.0, bool grid_entire_period = false,
            double n_sigma = 8.0, double kick_scale=1.0);
    /// Note: Use Space_charge_3d_open_hockney::get_internal_grid_shape for
    /// Distributed_fft3d.
    /// jfa: unnecessary complication
//    Space_charge_3d_open_hockney(Distributed_fft3d_sptr distributed_fft3d_sptr,
//            bool longitudinal_kicks = true, bool periodic_z = false,
//            double z_period = 0.0, bool grid_entire_period = false,
//            double n_sigma = 8.0);
    Space_charge_3d_open_hockney();
    virtual Space_charge_3d_open_hockney *
    clone();
    void
    set_diagnostics_list(Diagnostics_space_charge_3d_hockneys diagnosticss);
    void
    add_diagnostics(Diagnostics_space_charge_3d_hockney_sptr diagnostics_sptr);
    bool
    has_diagnostics();
    void
    setup_communication(Commxx_sptr const& bunch_comm_sptr);
    void
    setup_derived_communication();
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
    set_kick_scale(const double kick_scale);
    double
    get_kick_scale() const;

    void
    do_diagnostics(Rectangular_grid const& En, int component, double time_step, Step & step, Bunch & bunch);
    void
    apply_kick(Bunch & bunch, Rectangular_grid const& En, double delta_tau,
            int component);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger);
#endif
};

#endif /* SPACE_CHARGE_3D_OPEN_HOCKNEY_H_ */
