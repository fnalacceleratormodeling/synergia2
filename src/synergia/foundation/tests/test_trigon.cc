#include "synergia/utils/catch.hpp"

#include "synergia/foundation/trigon.h"

#include <complex>

template <unsigned int P> using Trig = Trigon<std::complex<double>, P, 1>;

TEST_CASE("factorials")
{
  static_assert(factorial(0) == 1);
  static_assert(factorial(1) == 1);
  static_assert(factorial(2) == 2);
  static_assert(factorial(3) == 6);
  static_assert(factorial(4) == 24);
  static_assert(factorial(5) == 120);
  static_assert(factorial(6) == 720);
  static_assert(factorial(7) == 5040);
  static_assert(factorial(8) == 40320);
  static_assert(factorial(9) == 362880);
  static_assert(factorial(10) == 3628800);
  static_assert(factorial(11) == 39916800);
  static_assert(factorial(12) == 479001600);
  // factorial(13) is too large to represent in a 32-bit integer, and so
  // it returns 0.
  static_assert(factorial(13) == 0);
}

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
