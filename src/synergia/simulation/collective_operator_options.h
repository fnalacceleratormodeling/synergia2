#ifndef COLLECTIVE_OPERATOR_OPTIONS_H
#define COLLECTIVE_OPERATOR_OPTIONS_H

#include <cereal/types/polymorphic.hpp>

class Collective_operator;

struct CO_options {
  virtual ~CO_options() = default;
  virtual CO_options* clone() const = 0;
  virtual Collective_operator* create_operator() const = 0;

  template <class Archive>
  void
  serialize(Archive& ar)
  {}
};

template <typename Derived, typename Operator>
struct CO_base_options : public CO_options {
  CO_options*
  clone() const override
  {
    return new Derived(static_cast<Derived const&>(*this));
  }

  Collective_operator*
  create_operator() const override
  {
    return new Operator(static_cast<Derived const&>(*this));
  }

  template <class Archive>
  void
  serialize(Archive& ar)
  {
    ar(cereal::base_class<CO_options>(this));
  }
};

#include <cereal/archives/json.hpp>

#endif
