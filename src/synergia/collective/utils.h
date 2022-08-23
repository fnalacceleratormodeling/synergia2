#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <cfenv>

#define IDXPRINT 0

KOKKOS_INLINE_FUNCTION
void
get_leftmost_indices_offset(double pos,
                            double left,
                            double inv_cell_size,
                            int& idx,
                            double& off)
{
  std::fesetround(FE_DOWNWARD);
  double scaled_location = (pos - left) * inv_cell_size - 0.5;
  idx = Kokkos::nearbyint(scaled_location);
  off = scaled_location - idx;

#if IDXPRINT == 1
  std::cout << "particle with pos : " << pos << ", left : " << left
            << ", inv_cell_size : " << inv_cell_size << ", idx : " << idx
            << ", off : " << off << std::endl;
#endif

  return;
}
