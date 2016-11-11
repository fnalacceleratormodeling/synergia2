#include "synergia/collective/deposit.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/core_diagnostics.h"

#include <iostream>

/// Deposit charge using Cloud-in-Cell (CIC) algorithm.
/// The indices on the rho array are in an unusual order: [z][y][x],
/// so that the FFTW routines can distribute along the z-axis.
/// The resulting charge density has units C/m^3.
void
deposit_charge_rectangular_zyx(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first)
{
    MArray3d_ref rho(rho_grid.get_grid_points());
    Const_MArray2d_ref parts(bunch.get_local_particles());
    if (zero_first) {
        for (unsigned int i = 0; i < rho.shape()[0]; ++i) {
            for (unsigned int j = 0; j < rho.shape()[1]; ++j) {
                for (unsigned int k = 0; k < rho.shape()[2]; ++k) {
                    rho[i][j][k] = 0.0;
                }
            }
        }
    }
    std::vector<double > h(rho_grid.get_domain().get_cell_size());
    double weight0 = (bunch.get_real_num() / bunch.get_total_num())
            * bunch.get_particle_charge() * pconstants::e
            / (h[0] * h[1] * h[2]);
    int ix, iy, iz;
    double offx, offy, offz;
    // jfa: This is probably a premature optimization. Two versions of the
    // deposit loop -- one for periodic, one for non-periodic.
    if (rho_grid.get_domain().is_periodic()) {
        for (int n = 0; n < bunch.get_local_num(); ++n) {
            // domain doesn't know about xyz->zyx transformation, so we
            // do it in the order of arguments here
            rho_grid.get_domain().get_leftmost_indices_offsets(
                    parts[n][4], parts[n][2], parts[n][0], iz, iy, ix, offz,
                    offy, offx);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int cellx = ix + i;
                        int celly = iy + j;
                        if ((cellx >0) && (cellx < int(rho.shape()[2]-1))
                                && (celly > 0)
                                && (celly < int(rho.shape()[1]-1))) {
                            int cellz = iz + k;
                            if (cellz >= 0) {
                                cellz = cellz % rho.shape()[0];
                            } else {
                                int period = rho.shape()[0];
                                cellz = period - 1 - ((-cellz - 1) % period);
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
            // domain doesn't know about xyz->zyx transformation, so we
            // do it in the order of arguments here
            rho_grid.get_domain().get_leftmost_indices_offsets(
                    parts[n][4], parts[n][2], parts[n][0], iz, iy, ix, offz,
                    offy, offx);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int cellx = ix + i;
                        int celly = iy + j;
                        int cellz = iz + k;                                               
//AM!  make sure the charge distribution is zero at the edge of the grid
// THIS toghether with zero electric field at the edge of the grid is essential for a conservative approximation                      
                        if ((cellx >0) && (cellx < int(rho.shape()[2]-1))
                                && (celly > 0)
                                && (celly < int(rho.shape()[1]-1))
                                && (cellz > 0)
                                && (cellz < int(rho.shape()[0]-1))) {                                                     
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

void
deposit_charge_rectangular_xyz(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first)
{
// the particles close (i.e at a distance smaller than cell_size/2) to the grid edges are not deposited
// they should aslo not be kicked by the electric field  
    
    std::vector<double > h(rho_grid.get_domain().get_cell_size());
    double weight0 = (bunch.get_real_num() / bunch.get_total_num())
            * bunch.get_particle_charge() * pconstants::e
            / (h[0] * h[1] * h[2]);
            
    
    MArray3d_ref rho(rho_grid.get_grid_points());
    Const_MArray2d_ref parts(bunch.get_local_particles());
    if (zero_first) {    
        for (unsigned int i = 0; i < rho.shape()[0]; ++i) {
            for (unsigned int j = 0; j < rho.shape()[1]; ++j) {
                for (unsigned int k = 0; k < rho.shape()[2]; ++k) {
                    rho[i][j][k] = 0.0;
                }
            }
        }
    } 

            
    int ix, iy, iz;
    double offx, offy, offz;
    if (rho_grid.get_domain().is_periodic()) {
        for (int n = 0; n < bunch.get_local_num(); ++n) {
            if (rho_grid.get_domain().get_leftmost_indices_offsets(
                        parts[n][0], parts[n][2], parts[n][4], ix, iy, iz, offx, offy, offz)){
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        int cellx = ix + i;
                        int celly = iy + j;
                        for (int k = 0; k < 2; ++k) {
                            if ((cellx >0) && (cellx < int(rho.shape()[0]-1))
                                && (celly > 0)
                                && (celly < int(rho.shape()[1]-1))) {
                                  int cellz = iz + k;
                                  int period = rho.shape()[2];
                                  cellz = (cellz % period)*(cellz >= 0)+(period - 1 - ((-cellz - 1) % period))*(cellz < 0);
                                  double weight = weight0 * (1 - i - (1 - 2 * i) * offx) * 
                                              (1 - j - (1 - 2 * j) * offy) *
                                              (1 - k - (1 - 2 * k) * offz); 
                                  rho[cellx][celly][cellz] += weight;                             
                           }
                       }
                    }
                }
            } 
        }         
    }
    else{
        for (int n = 0; n < bunch.get_local_num(); ++n) {
            if (rho_grid.get_domain().get_leftmost_indices_offsets(
                        parts[n][0], parts[n][2], parts[n][4], ix, iy, iz, offx, offy, offz)){
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        for (int k = 0; k < 2; ++k) {
                            int cellx = ix + i;
                            int celly = iy + j;
                            int cellz = iz + k;
                            if ((cellx >0) && (cellx < int(rho.shape()[0]-1))
                                && (celly > 0)
                                && (celly < int(rho.shape()[1]-1))
                                && (cellz > 0)
                                && (cellz < int(rho.shape()[2]-1))) { 
                                      double weight = weight0 * (1 - i - (1 - 2 * i) * offx) * 
                                                  (1 - j - (1 - 2 * j) * offy) *
                                                  (1 - k - (1 - 2 * k) * offz); 
                                      rho[cellx][celly][cellz] += weight;                       
                            }
                        }                   
                    }
                }
            } 
        }         
    }    
   //  rho_grid.set_normalization(total_charge_per_cell_vol*(h[0] * h[1] * h[2]));
       
}

void
deposit_charge_rectangular_2d(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first)
{
    MArray2dc_ref rho_2dc(rho_grid.get_grid_points_2dc());
    MArray1d_ref rho_1d(rho_grid.get_grid_points_1d());
    Const_MArray2d_ref parts(bunch.get_local_particles());
    std::vector<int > grid_shape(3);
    grid_shape[0] = rho_2dc.shape()[0];
    grid_shape[1] = rho_2dc.shape()[1];
    grid_shape[2] = rho_1d.shape()[0];
    if (zero_first) {
        for (int i = 0; i < grid_shape[0]; ++i) {           // x
            for (int j = 0; j < grid_shape[1]; ++j) {       // y
                rho_2dc[i][j] = 0.0;
            }
        }
        for (int k = 0; k < grid_shape[2]; ++k) {            // z
            rho_1d[k] = 0.0;
        }
    }
    std::vector<double > h(rho_grid.get_domain().get_cell_size());
    double weight0 = (bunch.get_real_num() / bunch.get_total_num())
            * bunch.get_particle_charge() * pconstants::e
            / (h[0] * h[1]); // * h[2]);
    int ix, iy, iz;
    double offx, offy, offz;
    if (grid_shape[2] == 1) {
        double mean(Core_diagnostics::calculate_z_mean(bunch));
        double std(Core_diagnostics::calculate_z_std(bunch, mean));
        rho_1d[0] = bunch.get_total_num() / (std::sqrt(12.0) * std);
    }
    for (int n = 0; n < bunch.get_local_num(); ++n) {
        // no xyz->zyx transformation
        rho_grid.get_domain().get_leftmost_indices_offsets(
                parts[n][0], parts[n][2], parts[n][4], ix, iy, iz, offx,
                offy, offz);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                int cellx = ix + i;
                int celly = iy + j;
                if ((cellx >= 0) && (cellx < grid_shape[0]) && (celly >= 0)
                        && (celly < grid_shape[1])) {
                    double weight = weight0 * (1 - i - (1 - 2 * i)
                            * offx) * (1 - j - (1 - 2 * j) * offy);
                    rho_2dc[cellx][celly] += weight;
                }
            }
        }
        for (int k = 0; k < 2; ++k) {
            int cellz = iz + k;
            if ((grid_shape[2] > 1) && (cellz >= 0) && (cellz < grid_shape[2])) {
                double weight = (1 - k - (1 - 2 * k) * offz) / h[2];
                rho_1d[cellz] += weight;
            }
        }
    }
}

void
deposit_charge_rectangular_2d(Rectangular_grid & rho_grid,
        Raw_MArray2d & particle_bin, Bunch const& bunch, bool zero_first)
{
    MArray2dc_ref rho_2dc(rho_grid.get_grid_points_2dc());
    MArray1d_ref rho_1d(rho_grid.get_grid_points_1d());
    Const_MArray2d_ref parts(bunch.get_local_particles());
    std::vector<int > grid_shape(3);
    grid_shape[0] = rho_2dc.shape()[0];
    grid_shape[1] = rho_2dc.shape()[1];
    grid_shape[2] = rho_1d.shape()[0];
    if (zero_first) {
        for (int i = 0; i < grid_shape[0]; ++i) {           // x
            for (int j = 0; j < grid_shape[1]; ++j) {       // y
                rho_2dc[i][j] = 0.0;
            }
        }
        for (int k = 0; k < grid_shape[2]; ++k) {            // z
            rho_1d[k] = 0.0;
        }
    }
    std::vector<double > h(rho_grid.get_domain().get_cell_size());
    double weight0 = (bunch.get_real_num() / bunch.get_total_num())
            * bunch.get_particle_charge() * pconstants::e
            / (h[0] * h[1]); // * h[2]);
    int ix, iy, iz;
    double offx, offy, offz;
    if (grid_shape[2] == 1) {
        double mean(Core_diagnostics::calculate_z_mean(bunch));
        double std(Core_diagnostics::calculate_z_std(bunch, mean));
        rho_1d[0] = bunch.get_local_num() / (std::sqrt(12.0) * std);
    }
    for (int n = 0; n < bunch.get_local_num(); ++n) {
        // no xyz->zyx transformation
        rho_grid.get_domain().get_leftmost_indices_offsets(
                parts[n][0], parts[n][2], parts[n][4], ix, iy, iz, offx,
                offy, offz);
        particle_bin.m[n][0] = ix;
        particle_bin.m[n][1] = offx;
        particle_bin.m[n][2] = iy;
        particle_bin.m[n][3] = offy;
        particle_bin.m[n][4] = iz;
        particle_bin.m[n][5] = offz;
        int cellx1, cellx2, celly1, celly2;
        cellx1 = ix;
        cellx2 = cellx1 + 1;
        celly1 = iy;
        celly2 = celly1 + 1;
        if ((cellx1 >= 0) && (cellx2 < grid_shape[0]) && (celly1 >= 0)
                && (celly2 < grid_shape[1])) {
            double aoffx, aoffy;
            aoffx = 1. - offx;
            aoffy = 1. - offy;
            rho_2dc[cellx1][celly1] += weight0 * aoffx * aoffy;
            rho_2dc[cellx1][celly2] += weight0 * aoffx * offy;
            rho_2dc[cellx2][celly1] += weight0 * offx * aoffy;
            rho_2dc[cellx2][celly2] += weight0 * offx * offy;
        }

        for (int k = 0; k < 2; ++k) {
            int cellz = iz + k;
            if ((grid_shape[2] > 1) && (cellz >= 0) && (cellz < grid_shape[2])) {
                double weight = (1 - k - (1 - 2 * k) * offz) / h[2];
                rho_1d[cellz] += weight;
            }
        }
    }
}
