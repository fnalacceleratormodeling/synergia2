#ifndef LINEAR_MATCHING_H_
#define LINEAR_MATCHING_H_

#include "synergia/bunch/bunch.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/distribution.h"

namespace linear_matching
{
    MArray2d
    get_covariance(Lattice_simulator & lattice_simulator, double hor_rms_emit,
            double vert_rms_emit, double long_rms_emit);

    void
    fill(Bunch & bunch, Lattice_simulator & lattice_simulator,
            double hor_rms_emit, double vert_rms_emit, double long_rms_emit,
            Distribution & distribution);

    void
    fill(Bunch & bunch, Lattice_simulator & lattice_simulator,
            double hor_rms_emit, double vert_rms_emit, double long_rms_emit,
            unsigned long int seed);

    void
    fill_transverse(Bunch & bunch, Lattice_simulator & lattice_simulator,
            double hor_rms_emit, double vert_rms_emit, double z_param,
            double dpop_rms, Distribution & distribution, bool uniform_z);

    void
    fill_transverse(Bunch & bunch, Lattice_simulator & lattice_simulator,
            double hor_rms_emit, double vert_rms_emit, double z_param,
            double dpop_rms, unsigned long int seed, bool uniform_z);
}

#endif /* LINEAR_MATCHING_H_ */
