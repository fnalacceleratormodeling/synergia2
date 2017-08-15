#include "linear_matching.h"

#include "Eigen/Eigen"
#include "Eigen/Cholesky"

using namespace Eigen;

MArray2d
linear_matching::get_covariance(Lattice_simulator & lattice_simulator,
        double hor_rms_emit, double vert_rms_emit, double long_rms_emit)
{

}

#if 0
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
#endif

