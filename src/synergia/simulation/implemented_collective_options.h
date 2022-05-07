#ifndef IMPLEMENTED_COLLECTIVE_OPTIONS_H
#define IMPLEMENTED_COLLECTIVE_OPTIONS_H

#include <array>
#include <memory>

#include "synergia/utils/cereal.h"

enum class green_fn_t {
  pointlike,
  linear,
};

enum class LongitudinalDistribution {
  gaussian,
  uniform,
};

struct Space_charge_3d_open_hockney_options {

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
    ar(shape);
    ar(doubled_shape);
    ar(green_fn);
    ar(periodic_z);
    ar(z_period);
    ar(grid_entire_period);
    ar(n_sigma);
    ar(comm_group_size);
  };
};

struct Space_charge_2d_open_hockney_options {

  std::array<int, 3> shape;
  std::array<int, 3> doubled_shape;
  bool periodic_z;
  double z_period;
  bool grid_entire_period;
  double n_sigma;
  int comm_group_size;

  Space_charge_2d_open_hockney_options(int gridx = 32,
                                       int gridy = 32,
                                       int gridz = 32)
    : shape{gridx, gridy, gridz}
    , doubled_shape{gridx * 2, gridy * 2, gridz}
    , periodic_z(false)
    , z_period(0.0)
    , grid_entire_period(false)
    , n_sigma(8.0)
    , comm_group_size(4)
  {}

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(shape);
    ar(doubled_shape);
    ar(periodic_z);
    ar(z_period);
    ar(grid_entire_period);
    ar(n_sigma);
    ar(comm_group_size);
  }
};

struct Space_charge_2d_kv_options {
  using LD = LongitudinalDistribution;

  // switch to control whether the linear charge density is assumed to be
  // distributed gaussian or uniform over the bunch length.
  LD longitudinal_distribution = LD::uniform;

  bool strictly_linear = true;
  bool strictly_centered = false;

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(longitudinal_distribution);
    ar(strictly_linear);
    ar(strictly_centered);
  }
};

struct Space_charge_rectangular_options {

  std::array<int, 3> shape;
  std::array<double, 3> pipe_size;
  int comm_group_size;

  Space_charge_rectangular_options(
    std::array<int, 3> const& shape = {32, 32, 64},
    std::array<double, 3> const& pipe_size = {0.1, 0.1, 1.0})
    : shape(shape), pipe_size(pipe_size), comm_group_size(1)
  {}

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(shape);
    ar(pipe_size);
    ar(comm_group_size);
  }
};

struct Impedance_options {

  std::string wake_file;
  std::string wake_type;
  int z_grid;
  bool full_machine;
  int nstored_turns;
  int num_buckets;
  double orbit_length;
  double bunch_spacing;
  std::array<int, 3> wn;

  Impedance_options(std::string const& wake_file = "",
                    std::string const& wake_type = "",
                    int z_grid = 1000)
    : wake_file(wake_file)
    , wake_type(wake_type)
    , z_grid(z_grid)
    , full_machine(false)
    , nstored_turns(15)
    , num_buckets(1)
    , orbit_length(1)
    , bunch_spacing(1)
  {}

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(wake_file);
    ar(wake_type);
    ar(z_grid);
    ar(full_machine);
    ar(nstored_turns);
    ar(num_buckets);
    ar(orbit_length);
    ar(bunch_spacing);
  }
};

struct Dummy_CO_options {

  template <class Archive>
  void
  serialize(Archive& ar)
  {}
};

#endif
