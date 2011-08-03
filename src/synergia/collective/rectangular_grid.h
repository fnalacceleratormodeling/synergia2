#ifndef RECTANGULAR_GRID_H_
#define RECTANGULAR_GRID_H_
#include <iostream>
#include <fstream>
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/utils/multi_array_typedefs.h"

class Rectangular_grid
{
private:
    Rectangular_grid_domain_sptr domain_sptr;
    boost::shared_ptr<MArray3d > grid_points_sptr;
    double normalization;
public:
    Rectangular_grid(std::vector<double > const & physical_size, std::vector<
            double > const & physical_offset,
            std::vector<int > const & grid_shape, bool periodic_z);
    Rectangular_grid(
            Rectangular_grid_domain_sptr rectangular_grid_domain_sptr);
    Rectangular_grid_domain_sptr
    get_domain_sptr() const;
    Rectangular_grid_domain_sptr
    get_domain_sptr();
    MArray3d_ref const&
    get_grid_points() const;
    MArray3d_ref &
    get_grid_points();
    void
    set_normalization(double val);
    double
    get_normalization() const;
    //
    // P.L. addition, Aug 3 2011 
    // 
    inline double get_interpolated(std::vector<double> location) const {
      return get_interpolated_coord(location[0], location[1], location[2]);
    } 
    //
    // Tested for the 2D case only, where coordinate X is X and Y is Y. iz = 0, for this case.
    // Other cases needs work most likely. But I am always confused on the convention in 
    // index name, Z vs X.  Depends on usage... 
    // 
    inline double get_interpolated_coord(double x, double y, double z) const {
       // tri-linear interpolation
       int ix, iy, iz;
       double offx, offy, offz;
       this->get_domain_sptr()->get_leftmost_indices_offsets(x, y, z, ix, iy, iz,
            offx, offy, offz);
       MArray3d_ref a(this->get_grid_points());
       double val;
       if ((get_domain_sptr()->get_grid_shape()[0] > 1) && (get_domain_sptr()->get_grid_shape()[1] > 1) &&
	      (get_domain_sptr()->get_grid_shape()[2] > 1)) {
       if ((ix < 0) || (ix >= get_domain_sptr()->get_grid_shape()[0] - 1) || (iy
              < 0) || (iy >= get_domain_sptr()->get_grid_shape()[1] - 1) || (iz
              < 0) || (iz >= get_domain_sptr()->get_grid_shape()[2] - 1)) {
          val = 0.0;
       } else {
          val = ((1.0 - offz) * (1.0 - offy) * (1.0 - offx) * a[iz][iy][ix]
                + offz * (1.0 - offy) * (1.0 - offx) * a[iz + 1][iy][ix] + (1.0
                - offz) * offy * (1.0 - offx) * a[iz][iy + 1][ix]
                + (1.0 - offz) * (1.0 - offy) * offx * a[iz][iy][ix + 1] + offz
                * offy * (1.0 - offx) * a[iz + 1][iy + 1][ix] + offz * (1.0
                - offy) * offx * a[iz + 1][iy][ix + 1] + (1.0 - offz) * offy
                * offx * a[iz][iy + 1][ix + 1] + offz * offy * offx
                * a[iz + 1][iy + 1][ix + 1]);
         }
       } else if (get_domain_sptr()->get_grid_shape()[0] == 1) {
	    // 2D,  Y-Z plane 
            if ((iy < 0) || (iy >= get_domain_sptr()->get_grid_shape()[1] - 1) || (iz
                    < 0) || (iz >= get_domain_sptr()->get_grid_shape()[2] - 1)) {
              val = 0.0;
            } else {
               val = ((1.0 - offz) * (1.0 - offy) * a[ix][iy][iz]
		 + offy * (1.0 - offz) * a[ix][iy + 1][iz]
                 + (1.0 - offy) * offz * a[ix][iy][iz + 1] 
		 + offy * offz * a[ix][iy + 1][iz + 1]);
            }
      } else if (get_domain_sptr()->get_grid_shape()[1] == 1) {
	    // 2D,  X-Z plane 
            if ((ix < 0) || (ix >= get_domain_sptr()->get_grid_shape()[0] - 1) || (iz
                    < 0) || (iz >= get_domain_sptr()->get_grid_shape()[2] - 1)) {
              val = 0.0;
            } else {
               val = ((1.0 - offz) * (1.0 - offx) * a[ix][iy][iz]
		 + offx * (1.0 - offz) * a[ix+1][iy][iz]
                 + (1.0 - offx) * offz * a[ix][iy][iz + 1] 
		 + offx * offz * a[ix + 1][iy][iz + 1]);
            } 
	    
       }  else if (get_domain_sptr()->get_grid_shape()[2] == 1) {
	    // 2D,  X-Y plane 
            if ((ix < 0) || (ix >= get_domain_sptr()->get_grid_shape()[0] - 1) || (iy
                    < 0) || (iy >= get_domain_sptr()->get_grid_shape()[1] - 1)) {
              val = 0.0;
            } else {
               val = ((1.0 - offy) * (1.0 - offx) * a[ix][iy][iz]
		 + offx * (1.0 - offy) * a[ix+1][iy][iz]
                 + (1.0 - offx) * offy * a[ix][iy+1][iz] 
		 + offx * offy * a[ix + 1][iy+1][iz]);
            }
       }
    return val;
     
    
    }
    
};

typedef boost::shared_ptr<Rectangular_grid> Rectangular_grid_sptr;

#endif /* RECTANGULAR_GRID_H_ */
