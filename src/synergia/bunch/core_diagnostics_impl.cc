#include "core_diagnostics.h"

karray1d
Core_diagnostics::calculate_mean(Bunch const& bunch)
{
  karray1d mean("mean", 6);
  const auto particles = bunch.get_local_particles();
  const auto masks = bunch.get_local_particle_masks();

  const auto total_bunch_particles = bunch.get_total_num();
  const auto local_bunch_capacity = bunch.size();

  auto instances = Kokkos::Experimental::partition_space(
    Kokkos::DefaultExecutionSpace(), 1, 1, 1, 1, 1, 1);

  return mean;
}
