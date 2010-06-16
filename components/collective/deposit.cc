#include "components/collective/deposit.h"

/// Deposit charge using Cloud-in-Cell (CIC) algorithm.
/// The indices on the rho array are in an unusual order: [z][y][x],
/// so that the FFTW routines can distribute along the z-axis.
void
deposit_charge_rectangular(Rectangular_grid & rho_grid, Bunch & bunch,
        bool zero_first)
{
    MArray3d_ref rho(rho_grid.get_grid_points());
    MArray2d_ref parts(bunch.get_local_particles());
    if (zero_first) {
        for (int i = 0; i < rho.shape()[0]; ++i) {
            for (int j = 0; j < rho.shape()[1]; ++j) {
                for (int k = 0; k < rho.shape()[2]; ++k) {
                    rho[i][j][k] = 0.0;
                }
            }
        }
    }
    std::vector<double > h(rho_grid.get_domain_sptr()->get_cell_size());
    double weight0 = 1.0 / (h[0] * h[1] * h[2]);
    int ix, iy, iz;
    double offx, offy, offz;
    // jfa: This is probably a premature optimization. Two versions of the
    // deposit loop -- one for periodic, one for non-periodic.
    if (rho_grid.get_domain_sptr()->is_periodic()) {
        std::cout << "jfa: deposit periodic\n";
        for (int n = 0; n < bunch.get_local_num(); ++n) {
            rho_grid.get_domain_sptr()->get_leftmost_indices_offsets(
                    parts[n][0], parts[n][2], parts[n][4], ix, iy, iz, offx,
                    offy, offz);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int cellx = ix + i;
                        int celly = iy + j;
                        if ((cellx >= 0) && (cellx < rho.shape()[2]) && (celly
                                >= 0) && (celly < rho.shape()[1])) {
                            int cellz = iz + k;
                            if (cellz >= 0) {
                                cellz = cellz % rho.shape()[0];
                            } else {
                                cellz = rho.shape()[0] - 1 + ((cellz + 1)
                                        % rho.shape()[0]);
                            }
                            double weight = weight0 * (1 - i - (1 - 2 * i)
                                    * offx) * (1 - j - (1 - 2 * j) * offy) * (1
                                    - k - (1 - 2 * k) * offz);
                            rho[cellz][celly][cellx] += weight;
                        }
                    }
                }
            }
        }
    } else {
        for (int n = 0; n < bunch.get_local_num(); ++n) {
            rho_grid.get_domain_sptr()->get_leftmost_indices_offsets(
                    parts[n][0], parts[n][2], parts[n][4], ix, iy, iz, offx,
                    offy, offz);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int cellx = ix + i;
                        int celly = iy + j;
                        int cellz = iz + k;
                        if ((cellx >= 0) && (cellx < rho.shape()[2]) && (celly
                                >= 0) && (celly < rho.shape()[1]) && (cellz
                                >= 0) && (cellz < rho.shape()[0])) {
                            double weight = weight0 * (1 - i - (1 - 2 * i)
                                    * offx) * (1 - j - (1 - 2 * j) * offy) * (1
                                    - k - (1 - 2 * k) * offz);
                            rho[cellz][celly][cellx] += weight;
                        }
                    }
                }
            }
        }
    }
}
