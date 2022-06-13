#ifndef CORE_DIAGNOSTICS_H_
#define CORE_DIAGNOSTICS_H_

#include "synergia/bunch/bunch.h"

struct Core_diagnostics {
  static karray1d calculate_mean(Bunch const& bunch);

  static double calculate_z_mean(Bunch const& bunch);

  static karray1d calculate_abs_mean(Bunch const& bunch);

  static karray1d calculate_std(Bunch const& bunch, karray1d const& mean);

  static karray2d_row calculate_sum2(Bunch const& bunch, karray1d const& mean);

  static karray2d_row calculate_mom2(Bunch const& bunch, karray1d const& mean);

  static karray1d calculate_min(Bunch const& bunch);

  static karray1d calculate_max(Bunch const& bunch);
};

#endif /* CORE_DIAGNOSTICS_H_ */
