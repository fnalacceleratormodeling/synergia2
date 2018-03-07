#include "synergia/collective/deposit.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/core_diagnostics.h"

#include "synergia/utils/synergia_omp.h"

/// Deposit charge using Cloud-in-Cell (CIC) algorithm.
/// The indices on the rho array are in an unusual order: [z][y][x],
/// so that the FFTW routines can distribute along the z-axis.
/// The resulting charge density has units C/m^3.
//
void
deposit_charge_rectangular_zyx(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first)
{
    MArray3d_ref rho(rho_grid.get_grid_points());
    Const_MArray2d_ref parts(bunch.get_local_particles());

    if (zero_first) {
        std::memset( rho.data(), 0, sizeof(double)*rho.num_elements() );
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


inline bool ingrid(int x, int gx)
{
  return x>=0 && x<gx;
}

inline bool ingrid(int x, int y, int z, int gx, int gy, int gz)
{
  return x>=0 && y>=0 && z>=0 && x<gx && y<gy && z<gz;
}
  

void
deposit_charge_rectangular_zyx_omp_reduce(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first)
{
    MArray3d_ref rho(rho_grid.get_grid_points());
    Const_MArray2d_ref parts(bunch.get_local_particles());

    const std::vector<double> & h    = rho_grid.get_domain_sptr()->get_cell_size();
    const std::vector<double> & offs = rho_grid.get_domain_sptr()->get_physical_offset();
    const std::vector<double> & size = rho_grid.get_domain_sptr()->get_physical_size();

    double w0 = (bunch.get_real_num() / bunch.get_total_num())
            * bunch.get_particle_charge() * pconstants::e
            / (h[0] * h[1] * h[2]);
 
    int npart = bunch.get_local_num();

    int gx = rho.shape()[2];
    int gy = rho.shape()[1];
    int gz = rho.shape()[0];

    int nc  = gx * gy * gz; // num of cells

    static int nt = 0;
    static int ncc = 0;
    static double * rl = 0;
    static double * rs = 0;

    if( nt==0 )
    {
        #pragma omp parallel
        { nt = omp_get_num_threads(); }

        ncc = nc;

        rl = new double[nt*nc];  // nt copies of +1 cells
        rs = new double[nc];  
    }

    if( nc!=ncc )
    {
        // bunch geometry has been changed
        delete [] rl;
        delete [] rs;

        ncc = nc;

        rl = new double[nt*nc];  // nt copies of +1 cells
        rs = new double[nc];  
    }

    double lx = offs[2] - size[2] / 2.0;
    double ly = offs[1] - size[1] / 2.0;
    double lz = offs[0] - size[0] / 2.0;

    double cx = h[2];
    double cy = h[1];
    double cz = h[0];

    // jfa: This is probably a premature optimization. Two versions of the
    // deposit loop -- one for periodic, one for non-periodic.
    if (rho_grid.get_domain().is_periodic()) {
#if 0
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
                        if ((cellx >= 0) && (cellx < int(rho.shape()[2]))
                                && (celly >= 0)
                                && (celly < int(rho.shape()[1]))) {
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
#endif
    } 
    else 
    {
        #pragma omp parallel shared(npart, parts, lx, ly, lz, cx, cy, cz, w0, gx, gy, gz, rl, rs, nc, rho)
        {
            int nt = omp_get_num_threads();
            int it = omp_get_thread_num();

            int np = npart;
            int le = npart / nt;
            int ps = it * le;
            int pe = (it == nt-1) ? np : (it+1)*le;

            // zero the worksheet
            std::memset(rl + it*nc, 0, sizeof(double)*nc);

            double scaled_location;
            double ox, oy, oz, w;
            int    ix, iy, iz;

            for(int n=ps; n<pe; ++n)
            {
                scaled_location = (parts[n][0] - lx) / cx - 0.5;
                ix = fast_int_floor(scaled_location);
                ox = scaled_location - ix;

                scaled_location = (parts[n][2] - ly) / cy - 0.5;
                iy = fast_int_floor(scaled_location);
                oy = scaled_location - iy;

                scaled_location = (parts[n][4] - lz) / cz - 0.5;
                iz = fast_int_floor(scaled_location);
                oz = scaled_location - iz;

                int base = iz * gx * gy;

                if( ingrid(ix  , iy  , iz, gx, gy, gz) ) 
                    rl[(base + (iy  )*gx + ix  )*nt+it] += w0*(1-ox)*(1-oy)*(1-oz);
                if( ingrid(ix+1, iy  , iz, gx, gy, gz) ) 
                    rl[(base + (iy  )*gx + ix+1)*nt+it] += w0*(  ox)*(1-oy)*(1-oz);
                if( ingrid(ix  , iy+1, iz, gx, gy, gz) ) 
                    rl[(base + (iy+1)*gx + ix  )*nt+it] += w0*(1-ox)*(  oy)*(1-oz);
                if( ingrid(ix+1, iy+1, iz, gx, gy, gz) ) 
                    rl[(base + (iy+1)*gx + ix+1)*nt+it] += w0*(  ox)*(  oy)*(1-oz);

                base = (iz+1) * gx * gy;

                if( ingrid(ix  , iy  , iz+1, gx, gy, gz) ) 
                    rl[(base + (iy  )*gx + ix  )*nt+it] += w0*(1-ox)*(1-oy)*(oz);
                if( ingrid(ix+1, iy  , iz+1, gx, gy, gz) ) 
                    rl[(base + (iy  )*gx + ix+1)*nt+it] += w0*(  ox)*(1-oy)*(oz);
                if( ingrid(ix  , iy+1, iz+1, gx, gy, gz) ) 
                    rl[(base + (iy+1)*gx + ix  )*nt+it] += w0*(1-ox)*(  oy)*(oz); 
                if( ingrid(ix+1, iy+1, iz+1, gx, gy, gz) ) 
                    rl[(base + (iy+1)*gx + ix+1)*nt+it] += w0*(  ox)*(  oy)*(oz); 
            }

            #pragma omp barrier

            // reduction
            le = gz/nt;
            ps = it*le;
            pe = (it == nt-1) ? gz : (ps + le);

            for(int z=ps; z<pe; ++z)
            {
                for(int y=0; y<gy; ++y)
                {
                    for(int x=0; x<gx; ++x)
                    {
                        w = 0.0;
                        for(int n=0; n<nt; ++n) w += rl[(z*gx*gy + y*gx + x)*nt+n];
                        rs[z*gx*gy + y*gx +x] = w;
                    }
                }
            }

            #pragma omp barrier

        } //  end of #pragma parallel

    } // end of periodic if

    memcpy( rho.data(), rs, sizeof(double)*nc );
}


void
deposit_charge_rectangular_zyx_omp_interleaved(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first)
{
    MArray3d_ref rho(rho_grid.get_grid_points());
    Const_MArray2d_ref parts(bunch.get_local_particles());

    int gx = rho.shape()[2];
    int gy = rho.shape()[1];
    int gz = rho.shape()[0];

    int npart = bunch.get_local_num();

    double * po = new double[npart*3];  // offx, offy, offz
    int    * pi = new int[npart];       // cell index of particles
    int    * pc = new int[gx*gy*gz+1];  // accumulated particle count in cells
    int    * count = new int[gx*gy*gz]; // particle count in cells
    int    * pll = new int[npart];      // indexed list

    static int nt = 0;

    if( nt==0 )
    {
        #pragma omp parallel
        { nt = omp_get_num_threads(); }
    }

    #pragma omp parallel for
    for(int i=0; i<gx*gy*gz+1; ++i) pc[i] = 0;
    //memset( pc   , 0, sizeof(int)*(gx*gy*gz+1) );

    #pragma omp parallel for
    for(int i=0; i<gx*gy*gz; ++i) count[i] = 0;
    //memset( count, 0, sizeof(int)*(gx*gy*gz) );

    if (zero_first) {
      #pragma omp parallel for
      for(int i=0; i<gx*gy*gz; ++i)
        ((double*)&rho[0][0][0])[i] = 0.0;
      //memset( &rho[0][0][0], 0, sizeof(double)*gx*gy*gz );
    }

    std::vector<double > h(rho_grid.get_domain_sptr()->get_cell_size());
    double weight0 = (bunch.get_real_num() / bunch.get_total_num())
            * bunch.get_particle_charge() * pconstants::e
            / (h[0] * h[1] * h[2]);
    int ix, iy, iz;
    double offx, offy, offz;

    const std::vector<double> & offs  = rho_grid.get_domain_sptr()->get_physical_offset();
    const std::vector<double> & size  = rho_grid.get_domain_sptr()->get_physical_size();

    double lx = offs[2] - size[2] / 2.0;
    double ly = offs[1] - size[1] / 2.0;
    double lz = offs[0] - size[0] / 2.0;

    double cx = h[2];
    double cy = h[1];
    double cz = h[0];

    double w0 = weight0;

    // jfa: This is probably a premature optimization. Two versions of the
    // deposit loop -- one for periodic, one for non-periodic.
    if (rho_grid.get_domain_sptr()->is_periodic()) 
    {
#if 0
        for (int n = 0; n < bunch.get_local_num(); ++n) 
        {
            // domain doesn't know about xyz->zyx transformation, so we
            // do it in the order of arguments here
            rho_grid.get_domain_sptr()->get_leftmost_indices_offsets(
                    parts[n][4], parts[n][2], parts[n][0], iz, iy, ix, offz,
                    offy, offx);

            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int cellx = ix + i;
                        int celly = iy + j;

                        if ((cellx >= 0) && (cellx < int(rho.shape()[2]))
                                && (celly >= 0)
                                && (celly < int(rho.shape()[1]))) {
                            int cellz = iz + k;

                            if (cellz >= 0) {
                                cellz = cellz % rho.shape()[0];
                            } else {
                                int period = rho.shape()[0];
                                cellz = period - 1 - ((-cellz - 1) % period);
                            }

                            double weight = weight0 * (1 - i - (1 - 2 * i) * offx) 
                                                    * (1 - j - (1 - 2 * j) * offy) 
                                                    * (1 - k - (1 - 2 * k) * offz);

                            rho[cellz][celly][cellx] += weight;
                        }

                    } // end of k
                } // end of j
            } // end of i
        }
#endif
    } 
    else 
    {
        // decl. for private variables
        int n, i, c_idx, x, y, z;
        int tid, nthreads, seg;
        double scaled_location;
        double w, ox, oy, oz;
      
        int * p;
        double * ws0;
        double * ws1;

        // count
        #pragma omp parallel for \
            shared( pi, po ) \
            private(n, c_idx, p, ix, iy, iz, offx, offy, offz, scaled_location)
        for(n=0; n<npart; ++n )
        {
            scaled_location = (parts[n][0] - lx) / cx - 0.5;
            ix = fast_int_floor(scaled_location);
            offx = scaled_location - ix;

            scaled_location = (parts[n][2] - ly) / cy - 0.5;
            iy = fast_int_floor(scaled_location);
            offy = scaled_location - iy;

            scaled_location = (parts[n][4] - lz) / cz - 0.5;
            iz = fast_int_floor(scaled_location);
            offz = scaled_location - iz;

            if( ix<0 || iy<0 || iz<0 || ix>=gx || iy>=gy || iz>=gz )
            {
              pi[n] = -1;
            }
            else
            {
              c_idx = iz*gx*gy + iy*gx + ix;
              pi[n] = c_idx;

              po[n*3+0] = offx;
              po[n*3+1] = offy;
              po[n*3+2] = offz;

              // gcc build-in atomic add
              __sync_fetch_and_add(pc+c_idx+1, 1);
              
              // another choice is to use the omp atomic
              //#pragma omp atomic
              //++pc[c_idx+1];
            }
        }

        // accumulate particle counts
        for(i=1; i<=gx*gy*gz; ++i)
          pc[i] += pc[i-1];

        int idx, pos;

        // build indexed list in parallel
        #pragma omp parallel for shared(pi, count) private(i, idx, pos)
        for(i=0; i<npart; ++i)
        {
          idx = pi[i];
          if( idx==-1 ) continue;
          pos = pc[idx] + __sync_fetch_and_add( count+idx, 1 );
          pll[pos] = i;
        }

        // deposit to cells
        #pragma omp parallel \
            shared( rho, pi, po, pc, count, pll  ) \
            private(x, y, z, ox, oy, oz, n, i, c_idx, w, ws0, ws1, tid, nthreads, seg) 
        {
          tid = omp_get_thread_num();
          nthreads = omp_get_num_threads();

          // temp work table to hold intermediate results
          ws0 = new double[(gx+1)*(gy+1)]();
          ws1 = new double[(gx+1)*(gy+1)]();

          // interleaved partitioning along z-axis to have 
          // evenly balanced loads for each thread
          for(z=tid; z<gz; z+=nthreads)
          {
            memset(ws0, 0, sizeof(double)*(gx+1)*(gy+1));
            memset(ws1, 0, sizeof(double)*(gx+1)*(gy+1));

            for(x=0; x<gx; ++x)
            {
              for(y=0; y<gy; ++y)
              {
                c_idx = z*gx*gy + y*gx + x;

                for(n=0; n<count[c_idx]; ++n)
                {
                  i = pll[pc[c_idx]+n]; // index of the particle

                  ox = po[i*3+0]; oy = po[i*3+1]; oz = po[i*3+2];

                  w = w0 * (1-ox) * (1-oy) * (1-oz); ws0[(y  )*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (1-oy) * (1-oz); ws0[(y  )*(gx+1) + x+1] += w;
                  w = w0 * (1-ox) * (  oy) * (1-oz); ws0[(y+1)*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (  oy) * (1-oz); ws0[(y+1)*(gx+1) + x+1] += w;
                
                  w = w0 * (1-ox) * (1-oy) * (  oz); ws1[(y  )*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (1-oy) * (  oz); ws1[(y  )*(gx+1) + x+1] += w;
                  w = w0 * (1-ox) * (  oy) * (  oz); ws1[(y+1)*(gx+1) + x  ] += w;
                  w = w0 * (  ox) * (  oy) * (  oz); ws1[(y+1)*(gx+1) + x+1] += w;
                } 

              } //end of y
            } // end of x
        
            // write to global memory
            for(x=0; x<gx; ++x)
              for(y=0; y<gy; ++y)
                rho[z][y][x] += ws0[ y*(gx+1) + x ];

            // a synchronization is needed to avoid data pollution
            #pragma omp barrier

            if( z!=gz-1 )
              for(x=0; x<gx; ++x)
                for(y=0; y<gy; ++y)
                  rho[z+1][y][x] += ws1[ y*(gx+1) + x ];

          } // end of z loop


          delete [] ws0;
          delete [] ws1;
        } 

    }

    delete [] po;
    delete [] pi;
    delete [] pc;
    delete [] count;
    delete [] pll;
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

void
deposit_charge_rectangular_2d_omp_reduce(Rectangular_grid & rho_grid,
        Raw_MArray2d & particle_bin, Bunch const& bunch, bool zero_first)
{
    MArray2dc_ref rho_2dc(rho_grid.get_grid_points_2dc());
    MArray1d_ref rho_1d(rho_grid.get_grid_points_1d());
    Const_MArray2d_ref parts(bunch.get_local_particles());

    int g0, g1, g2;
    g0 = rho_2dc.shape()[0];
    g1 = rho_2dc.shape()[1];
    g2 = rho_1d.shape()[0];

#if 0
    int G0 = g0 + 2;
    int G1 = g1 + 2;
    int G2 = g2 + 2;
#endif

    int npart = bunch.get_local_num();

    int nt = 1;

    #pragma omp parallel shared(nt)
    { nt = omp_get_num_threads(); }

    double * lrho2d = new double[g0*g1*nt];
    double * lrho1d = new double[g2*nt];

    if (zero_first) 
    {
        for (unsigned int i = 0; i < g0; ++i)         // x
        {  
            for (unsigned int j = 0; j < g1; ++j)     // y
            {
                rho_2dc[i][j] = 0.0;
            }
        }

        for (unsigned int k = 0; k < g2; ++k)         // z
        {
            rho_1d[k] = 0.0;
        }
    }

    std::vector<double > h(rho_grid.get_domain().get_cell_size());
    double weight0 = (bunch.get_real_num() / bunch.get_total_num())
            * bunch.get_particle_charge() * pconstants::e
            / (h[0] * h[1]); // * h[2]);

    if (g2 == 1) 
    {
        double mean(Core_diagnostics::calculate_z_mean(bunch));
        double std(Core_diagnostics::calculate_z_std(bunch, mean));
        rho_1d[0] = bunch.get_local_num() / (std::sqrt(12.0) * std);
    }

    #pragma omp parallel \
        default(none) \
        shared(nt, npart, parts, lrho2d, lrho1d, h, \
               g0, g1, g2, particle_bin, \
               weight0, rho_2dc, rho_1d, rho_grid)
    {
        int ix, iy, iz;
        double offx, offy, offz;

        int it = omp_get_thread_num();

        int l  = npart / nt;                  // length
        int s  = it * l;                      // start particle
        int e  = (it==nt-1) ? npart : (s+l);  // end particle

        double * r2d = lrho2d + it * g0 * g1;
        double * r1d = lrho1d + it * g2;

        // zero buffer
        std::memset( r2d, 0, sizeof(double)*g0*g1 );
        std::memset( r1d, 0, sizeof(double)*g2    );

        for (int n = s; n < e; ++n) {
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

            int cellz1 = iz;
            int cellz2 = cellz1 + 1;

            if( cellz1>=0 && cellz1<g2 ) r1d[cellz1] += (1.0 - offz) / h[2];
            if( cellz2>=0 && cellz2<g2 ) r1d[cellz2] += offz / h[2];

            if( ix<0 || ix>g0-1 || iy<0 || iy>g1-1 ) continue;

            int cellx1, cellx2, celly1, celly2;
            cellx1 = ix;
            cellx2 = ix + 1;
            celly1 = iy;
            celly2 = iy + 1;

            double aoffx, aoffy;
            aoffx = 1. - offx;
            aoffy = 1. - offy;

            r2d[celly1*g0 + cellx1] += weight0 * aoffx * aoffy;
            r2d[celly1*g0 + cellx2] += weight0 *  offx * aoffy;
            r2d[celly2*g0 + cellx1] += weight0 * aoffx *  offy;
            r2d[celly2*g0 + cellx2] += weight0 *  offx *  offy;
        }

        // set boundary to zero
        for (int x=0; x<g0; ++x) 
        {
            r2d[x] = 0.0;
            r2d[(g1-1)*g0 + x] = 0.0;
        }

        for (int y=0; y<g1; ++y)
        {
            r2d[y*g0] = 0.0;
            r2d[y*g0 + g0 - 1] = 0.0;
        }
    }

    for (int t = 0; t < nt; ++t)
    {
        for (int y = 0; y < g1; ++y)
            for (int x = 0; x < g0; ++x)
                rho_2dc[x][y] += lrho2d[t*g0*g1 + y*g0 + x];

        if (g2 > 1)
            for (int z = 0; z < g2; ++z)
                rho_1d[z] += lrho1d[t*g2 + z];
    }

    delete [] lrho2d;
    delete [] lrho1d;
}
