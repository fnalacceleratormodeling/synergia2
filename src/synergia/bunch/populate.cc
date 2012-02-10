#include "populate.h"
#include "diagnostics.h"

#include "synergia/utils/eigen2/Eigen/Eigen"
#include "synergia/utils/eigen2/Eigen/Cholesky"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

#include "synergia/utils/multi_array_print.h"

void
adjust_moments(Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Core_diagnostics::calculate_mom2(bunch, bunch_mean));
    Matrix<double, 6, 6, Eigen::RowMajor > C(covariances.origin());
    Matrix<double, 6, 6, Eigen::RowMajor > X(bunch_mom2.origin());
    Matrix<double, 6, 6, Eigen::RowMajor > A = C.llt().matrixL()
            * X.llt().matrixL().inverse();

    int num_particles = bunch.get_local_num();
    Eigen::Map<Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >
            rho7(bunch.get_local_particles().origin(), num_particles, 7);
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

void
populate_6d(Distribution &dist, Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    MArray2d_ref particles(bunch.get_local_particles());
    for (int i = 0; i < 6; ++i) {
        dist.fill_unit_gaussian(particles[boost::indices[range()][i]]);
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
populate_transverse_KV_GaussLong(Distribution &dist, Bunch &bunch, double epsilMax,
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
      const double a2X = std::sqrt((x1*x1 + x2*x2) * epsilMax*beta_x); // The amplitude.. X physical plane
      particles[part][Bunch::x] = a2X*std::sin(phi2X);
      particles[part][Bunch::xp] = (1.0/beta_x)*(a2X*std::cos(phi2X) - alpha_x*particles[part][Bunch::x]);
      //  Repeat in y,y' plane 
      const double phi2Y = std::atan2(w,z);
      const double a2Y = std::sqrt((w*w + z*z)*epsilMax*beta_y); // The amplitude.. 
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
