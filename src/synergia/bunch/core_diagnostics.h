#ifndef CORE_DIAGNOSTICS_H_
#define CORE_DIAGNOSTICS_H_

#include "synergia/bunch/bunch.h"

struct Core_diagnostics
{
    static MArray1d
    calculate_mean(Bunch const& bunch);

    static double
    calculate_z_mean(Bunch const& bunch);

    static MArray1d
    calculate_spatial_mean(Bunch const& bunch);

    static MArray1d
    calculate_std(Bunch const& bunch, MArray1d_ref const& mean);

    static MArray1d
    calculate_spatial_std(Bunch const& bunch, MArray1d_ref const& mean);

    static MArray2d
    calculate_mom2(Bunch const& bunch, MArray1d_ref const& mean);

    static MArray1d
    calculate_min(Bunch const& bunch);

    static MArray1d
    calculate_max(Bunch const& bunch);
};


#endif /* CORE_DIAGNOSTICS_H_ */
