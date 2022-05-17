#ifndef COLLECTIVE_OPERATOR_OPTIONS_H
#define COLLECTIVE_OPERATOR_OPTIONS_H

#include "synergia/simulation/collective_operator.h"
#include "synergia/simulation/dummy_collective_operator.h"
#include "synergia/simulation/implemented_collective_operators.h"
#include "synergia/simulation/implemented_collective_options.h"

#include <variant>

using CO_options = std::variant<Dummy_CO_options,
                                Space_charge_3d_open_hockney_options,
                                Space_charge_2d_open_hockney_options,
                                Space_charge_2d_kv_options,
                                Space_charge_rectangular_options,
                                Impedance_options>;

struct create_collective_operator {

  std::shared_ptr<Operator>
  operator()(const Dummy_CO_options& ops)
  {
    return std::make_shared<Dummy_collective_operator>();
  }

  std::shared_ptr<Operator>
  operator()(const Space_charge_3d_open_hockney_options& ops)
  {
    return std::make_shared<Space_charge_3d_open_hockney>(
      static_cast<const Space_charge_3d_open_hockney_options&>(ops));
  }

  std::shared_ptr<Operator>
  operator()(const Space_charge_2d_open_hockney_options& ops)
  {
    return std::make_shared<Space_charge_2d_open_hockney>(
      static_cast<const Space_charge_2d_open_hockney_options&>(ops));
  }

  std::shared_ptr<Operator>
  operator()(const Space_charge_2d_kv_options& ops)
  {
    return std::make_shared<Space_charge_2d_kv>(
      static_cast<const Space_charge_2d_kv_options&>(ops));
  }

  std::shared_ptr<Operator>
  operator()(const Space_charge_rectangular_options& ops)
  {
    return std::make_shared<Space_charge_rectangular>(
      static_cast<const Space_charge_rectangular_options&>(ops));
  }

  std::shared_ptr<Operator>
  operator()(const Impedance_options& ops)
  {
    return std::make_shared<Impedance>(
      static_cast<const Impedance_options&>(ops));
  }
};

#endif
