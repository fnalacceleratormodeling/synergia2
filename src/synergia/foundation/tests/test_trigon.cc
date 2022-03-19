#include "synergia/utils/catch.hpp"

#include "synergia/foundation/trigon.h"

#include <complex>

template <unsigned int P> using Trig = Trigon<std::complex<double>, P, 1>;

TEST_CASE("Power 0")
{
  using trig_t = Trig<0>;
  static_assert(std::is_same_v<trig_t::data_type, std::complex<double>>);
}

TEST_CASE("Power 1")
{
  using trig_t = Trig<1>;
  static_assert(std::is_same_v<trig_t::data_type, std::complex<double>>);
}
