#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>

#include "synergia/lattice/lattice_element_slice.h"

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_operation.h"

#include "synergia/utils/logger.h"

class Propagator;

class Operator {

private:
  std::string name;
  std::string type;

  double time_fraction;

private:
  virtual void create_operations_impl(Lattice const& lattice) = 0;

  virtual void apply_impl(Bunch_simulator& simulator,
                          double time_step,
                          Logger& logger) = 0;

  virtual void print_impl(Logger& logger) const = 0;

public:
  Operator(std::string const& name, std::string const& type, double time)
    : name(name), type(type), time_fraction(time)
  {}

  virtual ~Operator() = default;

  std::string const&
  get_name() const
  {
    return name;
  }
  std::string const&
  get_type() const
  {
    return type;
  }
  double
  get_time_fraction() const
  {
    return time_fraction;
  }

  void
  create_operations(Lattice const& lattice)
  {
    create_operations_impl(lattice);
  }

  void
  apply(Bunch_simulator& simulator, double time_step, Logger& logger)
  {
    apply_impl(simulator, time_step * time_fraction, logger);
  }

  void
  print(Logger& logger) const
  {
    logger(LoggerV::DEBUG) << "operator name = " << name << ", type = " << type
                           << ", time_fraction = " << time_fraction << "\n";

    print_impl(logger);
  }
};

// collective
class Collective_operator : public Operator {
private:
  void
  create_operations_impl(Lattice const& lattice) final
  {}
  void
  print_impl(Logger& logger) const override
  {}

public:
  Collective_operator(std::string const& name, double time)
    : Operator(name, "collective", time)
  {}
};

// Independent
class Independent_operator : public Operator {
private:
  void apply_impl(Bunch_simulator& simulator,
                  double time_step,
                  Logger& logger) override;

  void create_operations_impl(Lattice const& lattice) override;

  void print_impl(Logger& logger) const override;

  bool need_update(Reference_particle const& ref, Logger& logger);

  void update_operations(Reference_particle const& ref);

private:
  std::vector<Lattice_element_slice> slices;
  std::vector<std::unique_ptr<Independent_operation>> operations;

public:
  Independent_operator(std::string const& name, double time);

  template <class... Args>
  Independent_operator&
  append_slice(Args&&... args)
  {
    slices.emplace_back(std::forward<Args>(args)...);
    return *this;
  }

  std::vector<Lattice_element_slice> const&
  get_slices() const
  {
    return slices;
  }

  friend class Propagator;
};

#endif /* OPERATOR_H_ */
