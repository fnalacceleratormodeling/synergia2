#include "synergia/utils/catch.hpp"

#include "synergia/foundation/physical_constants.h"

#include "synergia/collective/deposit.h"

const double mass = 100.0;
const double total_energy = 125.0;

const int total_num = 100;
const double real_num = 2.0e12;
auto const id = Bunch::id;

TEST_CASE("NoParticles", "[NoParticles]")
{
  Four_momentum fm(mass, total_energy);
  Reference_particle ref(pconstants::proton_charge, fm);

  Bunch bunch(ref, 0, 0, Commxx());
  CHECK(bunch.get_local_num() == 0);

  Rectangular_grid_domain domain(
    {4, 4, 4}, {1e-6, 1e-6, 1e-6}, {0, 0, 0}, false);

  const std::array<int, 3> dims{8, 8, 8};
  karray1d_dev rho_dev("rho_dev", 512);

  deposit_charge_rectangular_3d_kokkos_scatter_view(
    rho_dev, domain, dims, bunch);

  karray1d_hst rho_dev_hst = Kokkos::create_mirror_view(rho_dev);
  Kokkos::deep_copy(rho_dev_hst, rho_dev);
  Kokkos::fence();

  for (int idx = 0; idx < 512; idx++) {
    CHECK(rho_dev_hst(0) == Approx(0).margin(.01));
  }
}

TEST_CASE("OneParticle", "[OneParticle]")
{
  Four_momentum fm(mass, total_energy);
  Reference_particle ref(pconstants::proton_charge, fm);

  Bunch bunch(ref, 1, 1, Commxx());
  CHECK(bunch.get_local_num() == 1);

  Rectangular_grid_domain domain(
    {3, 3, 3}, {1e-6, 1e-6, 1e-6}, {0, 0, 0}, false);

  const std::array<int, 3> dims{3, 3, 3};
  karray1d_dev rho_dev("rho_dev", 3 * 3 * 3);

  deposit_charge_rectangular_3d_kokkos_scatter_view(
    rho_dev, domain, dims, bunch);

  karray1d_hst rho_dev_hst = Kokkos::create_mirror_view(rho_dev);
  Kokkos::deep_copy(rho_dev_hst, rho_dev);
  Kokkos::fence();

  std::array<int, 8> deposit_locs{12, 13, 15, 16, 22, 23, 25, 26};

  auto sum = 0.0i;
  for (int idx = 0; idx < 27; idx++) { sum += rho_dev_hst(idx); }
  CHECK(sum == Approx(1).margin(.01));

  for (auto idx : deposit_locs) {
    CHECK(rho_dev_hst(idx) == Approx(0.125).margin(.01));
  }
}
