#include <sstream>
#include <stdexcept>
#include "populate.h"
#include "diagnostics.h"

#include "eigen3/Eigen/Eigen"
#include "eigen3/Eigen/Cholesky"

using namespace Eigen;

#include "synergia/utils/floating_point.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/multi_array_assert.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;

namespace {
bool is_symmetric66(Const_MArray2d_ref &m) {
  bool symmetric = true;
  const double tolerance = 1.0e-14;
  for (int i = 0; i < 6; ++i) {
    for (int j = i + 1; j < 6; ++j) {
      if (!floating_point_equal(m[i][j], m[j][i], tolerance)) {
        symmetric = false;
      }
    }
  }
  return symmetric;
}
}

void
adjust_moments(Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    if (!is_symmetric66(covariances)) {
        throw std::runtime_error("adjust_moments: covariance matrix must be symmetric");
    }
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Core_diagnostics::calculate_mom2(bunch, bunch_mean));
    Matrix<double, 6, 6, Eigen::ColMajor > C(covariances.origin());
    Matrix<double, 6, 6, Eigen::ColMajor > G(C.llt().matrixL());
    Matrix<double, 6, 6, Eigen::ColMajor > X(bunch_mom2.origin());
    Matrix<double, 6, 6, Eigen::ColMajor > H(X.llt().matrixL());
    Matrix<double, 6, 6, Eigen::ColMajor > A(G * H.inverse());
    // jfa: dummy exists only to work around a bad interaction betwen
    //      Eigen3 and g++ 4.1.2
    std::stringstream dummy;
    dummy << C;

    int num_particles = bunch.get_local_num();
    int num_particles_padded = bunch.get_local_num_padded();

    Eigen::Map<Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor > >
            rho7(bunch.get_local_particles().origin(), num_particles_padded, 7);
    Matrix<double, 1, 6 > rhobar6(bunch_mean.origin());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        rho7.block<1, 6 > (part, 0) -= rhobar6;
    }

    rho7.block(0, 0, num_particles, 6) *= A.transpose();

    Matrix<double, 1, 6 > means6(means.origin());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        rho7.block<1, 6 > (part, 0) += means6;
    }
}

namespace
{
    void
    fill_unit_6d(Distribution & dist, MArray2d_ref & particles,
            Const_MArray2d_ref covariances, int start, int end)
    {
        for (int j = 0; j < 6; ++j) {
            dist.fill_unit_gaussian(
                    particles[boost::indices[range(start, end)][j]]);
            double scale = sqrt(covariances[j][j]);
            for (int i = start; i < end; ++i) {
                particles[i][j] *= scale;
            }
        }
    }

    inline
    bool
    good(MArray2d_ref particles, Const_MArray1d_ref limits, int index)
    {
        bool retval = true;
        for (int i = 0; i < 6; ++i) {
            double val = particles[index][i];
            double limit = limits[i];
            if ((limit > 0) && ((val > limit) or (val < -limit))) {
                retval = false;
            }
        }
        return retval;
    }

    void
    strip_unit_6d(Bunch & bunch, Const_MArray1d_ref limits, int & total_num,
            int & local_num)
    {
        MArray2d_ref particles(bunch.get_local_particles());
        local_num = bunch.get_local_num();
        int index = 0;
        while (index < local_num) {
            if (good(particles, limits, index)) {
                ++index;
            } else {
                int last = local_num - 1;
                if (good(particles, limits, last)) {
                    particles[index][0] = particles[last][0];
                    particles[index][1] = particles[last][1];
                    particles[index][2] = particles[last][2];
                    particles[index][3] = particles[last][3];
                    particles[index][4] = particles[last][4];
                    particles[index][5] = particles[last][5];
                }
                --local_num;
            }
        }
        MPI_Allreduce(&local_num, &total_num, 1, MPI_INT, MPI_SUM,
                bunch.get_comm().get());
    }
}

void
populate_6d(Distribution &dist, Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    MArray1d limits(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        limits[i] = 0.0;
    }
    populate_6d_truncated(dist, bunch, means, covariances, limits);
}

void
populate_6d_truncated(Distribution &dist, Bunch &bunch,
        Const_MArray1d_ref means, Const_MArray2d_ref covariances,
        Const_MArray1d_ref limits)
{
    multi_array_assert_size(means, 6, "populate_6d: means");
    multi_array_assert_size(covariances, 6, 6, "populate_6d: covariances");
    multi_array_assert_size(limits, 6, "populate_6d: limits");
    MArray2d_ref particles(bunch.get_local_particles());
    MArray2d unit_covariances(boost::extents[6][6]);
    MArray1d zero_means(boost::extents[6]);
    bool truncated(false);
    for (int i = 0; i < 6; ++i) {
        double n = limits[i];
        if (n > 0) {
            truncated = true;
            double cutoff_integral = exp(-n * n / 2.0)
                    * (sqrt(pi) * exp(n * n / 2.0) * erf(n / sqrt(2.0))
                            - sqrt(2.0) * n) / (sqrt(pi));
            unit_covariances[i][i] = 1 / (cutoff_integral * cutoff_integral);
        } else {
            unit_covariances[i][i] = 1.0;
        }
    }
    int start = 0;
    int end = bunch.get_local_num();
    fill_unit_6d(dist, particles, unit_covariances, start, end);
    if (truncated) {
        adjust_moments(bunch, zero_means, unit_covariances);
        int iteration = 0;
        int total_num, local_num;
        strip_unit_6d(bunch, limits, total_num, local_num);
        while (total_num < bunch.get_total_num()) {
            ++iteration;
            const int max_iterations = 50;
            if (iteration > max_iterations) {
                throw std::runtime_error(
                        "populate_6d_truncated: maximum number of truncation iterations exceeded. Algorithm known to fail ~< 2.5 sigma.");
            }
            fill_unit_6d(dist, particles, unit_covariances, local_num, end);
            adjust_moments(bunch, zero_means, unit_covariances);
            strip_unit_6d(bunch, limits, total_num, local_num);
        }
    }
    adjust_moments(bunch, means, covariances);
    bunch.check_pz2_positive();
}

void
populate_transverse_gaussian(Distribution &dist, Bunch &bunch,
        Const_MArray1d_ref means, Const_MArray2d_ref covariances, double cdt)
{
    MArray2d_ref particles(bunch.get_local_particles());
    for (int i = 0; i < 4; ++i) {
        dist.fill_unit_gaussian(particles[boost::indices[range()][i]]);
    }
    dist.fill_uniform(particles[boost::indices[range()][4]], 0.0, 1.0);
    dist.fill_unit_gaussian(particles[boost::indices[range()][5]]);

    MArray1d means_modified(means);
    means_modified[4] = 0.0;

    // Symmetry requires no correlations with the cdt coordinate. Make a copy
    // of the covariance matrix and manually set all correlations to zero.
    MArray2d covariances_modified(covariances);
    for (int k = 0; k < 6; ++k) {
        covariances_modified[k][4] = covariances_modified[4][k] = 0.0;
    }
    covariances_modified[4][4] = cdt * cdt / 12.0;
    adjust_moments(bunch, means_modified, covariances_modified);
}

void
populate_uniform_cylinder(Distribution &dist, Bunch &bunch, double radius,
        double cdt, double stdxp, double stdyp, double stddpop)
{
    MArray2d_ref particles(bunch.get_local_particles());
    dist.fill_unit_disk(particles[boost::indices[range()][Bunch::x]],
            particles[boost::indices[range()][Bunch::y]]);
    dist.fill_uniform(particles[boost::indices[range()][Bunch::cdt]], -cdt
            / 2.0, cdt / 2.0);
    dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::xp]]);
    dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::yp]]);
    dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::dpop]]);

    for (int part = 0; part < bunch.get_local_num(); ++part) {
        particles[part][Bunch::x] *= radius;
        particles[part][Bunch::y] *= radius;
        particles[part][Bunch::xp] *= stdxp;
        particles[part][Bunch::yp] *= stdyp;
        particles[part][Bunch::dpop] *= stddpop;
    }
}
void
populate_transverse_KV_GaussLong(Distribution &dist, Bunch &bunch, double epsilMax_x, double epsilMax_y,
        double alpha_x, double beta_x, double alpha_y, double beta_y,
        double cdt, double stddpop){
    MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
//
// First part: generate uniformally distributed on a hypershere, dim=4
// Algorith by marsaglia, found in the book by Wolfram,
// see  http://mathworld.wolfram.com/HyperspherePointPicking.html
//
      double d12 = 2.;
      double x1 = 0.; double x2=0.;
// Pick first wo points within a 2D circle, flat distribution
      while (d12 > 1.0) {
        x1 = 2.0*dist.get() - 1.0;
        x2 = 2.0*dist.get() - 1.0; 
        d12 = std::sqrt(x1*x1 + x2*x2);
      }
// Two more       
      double d34 = 2.;
      double x3 = 0.; double x4=0.;
      while (d34 > 1.0) {
        x3 = 2.0*dist.get() - 1.0;
        x4 = 2.0*dist.get() - 1.0; 
        d34 = std::sqrt(x3*x3 + x4*x4);
      }
// The 4 points on the 4-sphere are x1, x2, z, w 
      const double z = x3 * std::sqrt((1.0 - x1*x1 - x2*x2)/(x3*x3 + x4*x4));
      const double w = x4 * std::sqrt((1.0 - x1*x1 - x2*x2)/(x3*x3 + x4*x4));
// Now move from normal coordinate to physical using lattice functios
      const double phi2X = std::atan2(x2,x1);
      const double a2X = std::sqrt((x1*x1 + x2*x2) * epsilMax_x*beta_x); // The amplitude.. X physical plane
      particles[part][Bunch::x] = a2X*std::sin(phi2X);
      particles[part][Bunch::xp] = (1.0/beta_x)*(a2X*std::cos(phi2X) - alpha_x*particles[part][Bunch::x]);
      //  Repeat in y,y' plane 
      const double phi2Y = std::atan2(w,z);
      const double a2Y = std::sqrt((w*w + z*z)*epsilMax_y*beta_y); // The amplitude.. 
      particles[part][Bunch::y] = a2Y*std::sin(phi2Y);
      particles[part][Bunch::yp] = (1.0/beta_y)*(a2Y*std::cos(phi2Y) - alpha_y*particles[part][Bunch::y]);
    }
    dist.fill_uniform(particles[boost::indices[range()][Bunch::cdt]], -cdt
            / 2.0, cdt / 2.0);
    dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::dpop]]);
    for (int part = 0; part < bunch.get_local_num(); ++part) particles[part][Bunch::dpop] *= stddpop;
}
void
populate_two_particles(Bunch &bunch,
         double p1x, double p1xp, double p1y, double p1yp, double p1cdt, double p1dpop, 
         double p2x, double p2xp, double p2y, double p2yp, double p2cdt, double p2dpop) {
    MArray2d_ref particles(bunch.get_local_particles());
    if (bunch.get_local_num() !=2) {
        std::ostringstream errMsgStream; errMsgStream << "Expecting only two particles when" 
	                                              << bunch.get_local_num() << "generated";
        std::string errMsg(errMsgStream.str());
        throw std::runtime_error(errMsg.c_str());
    }     
    particles[0][Bunch::x] = p1x; particles[0][Bunch::xp] = p1xp;
    particles[0][Bunch::y] = p1y; particles[0][Bunch::yp] = p1yp;
    particles[0][Bunch::cdt] = p1cdt; particles[0][Bunch::dpop] = p1dpop;
    particles[1][Bunch::x] = p2x; particles[1][Bunch::xp] = p2xp;
    particles[1][Bunch::y] = p2y; particles[1][Bunch::yp] = p2yp;
    particles[1][Bunch::cdt] = p2cdt; particles[1][Bunch::dpop] = p2dpop;
} 
