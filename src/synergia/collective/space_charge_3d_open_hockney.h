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
};

#endif /* SPACE_CHARGE_3D_OPEN_HOCKNEY_H_ */
