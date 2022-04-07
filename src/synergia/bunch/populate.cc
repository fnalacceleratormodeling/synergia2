#include <sstream>
#include <stdexcept>

#include "core_diagnostics.h"
#include "populate.h"
#include "populate_host.h"

#include "synergia/foundation/math_constants.h"
#include "synergia/utils/floating_point.h"

using mconstants::pi;

namespace {
  bool
  is_symmetric66(const_karray2d_row m)
  {
    bool symmetric = true;
    const double tolerance = 1.0e-14;

    for (int i = 0; i < 6; ++i) {
      for (int j = i + 1; j < 6; ++j) {
        if (!floating_point_equal(m(i, j), m(j, i), tolerance)) {
          symmetric = false;
        }
      }
    }

    return symmetric;
  }
}

void
adjust_moments(Bunch& bunch,
               const_karray1d means,
               const_karray2d_row covariances)
{
  if (!is_symmetric66(covariances))
    throw std::runtime_error(
      "adjust_moments: covariance matrix must be symmetric");

  // calculate_mean and mom2 are performed on device memory, so we need to
  // copy the particle data from host to device first. checkout is not
  // necessary since the core diagnostics do not change the particle data
  bunch.checkin_particles();

  karray1d bunch_mean = Core_diagnostics::calculate_mean(bunch);
  karray2d_row bunch_mom2 = Core_diagnostics::calculate_mom2(bunch, bunch_mean);

  int num_particles = bunch.size();
  int num_particles_slots = bunch.capacity();

  // auto strides = bunch.get_particle_strides();
  // int num_particles_slots = strides[1];

  adjust_moments_host(means.data(),
                      covariances.data(),
                      bunch_mean.data(),
                      bunch_mom2.data(),
                      num_particles,
                      num_particles_slots,
                      bunch.get_host_particles().data());
}

namespace {
  void
  fill_unit_6d(Distribution& dist,
               HostParticles particles,
               const_karray2d_row covariances,
               int start,
               int end)
  {
    for (int j = 0; j < 6; ++j) {
      const double scale = sqrt(covariances(j, j));

      for (int p = start; p < end; ++p)
        particles(p, j) = dist.get_unit_gaussian() * scale;
    }
  }

  inline bool
  good(ConstHostParticles particles, const_karray1d limits, int index)
  {
    bool retval = true;

    for (int i = 0; i < 6; ++i) {
      double val = particles(index, i);
      double limit = limits[i];

      if ((limit > 0) && ((val > limit) or (val < -limit))) retval = false;
    }

    return retval;
  }

  void
  strip_unit_6d(Bunch& bunch,
                const_karray1d limits,
                int& total_num,
                int& local_num)
  {
    auto particles = bunch.get_host_particles();
    local_num = bunch.get_local_num();

    int index = 0;
    while (index < local_num) {
      if (good(particles, limits, index)) {
        ++index;
      } else {
        int last = local_num - 1;
        if (good(particles, limits, last)) {
          particles(index, 0) = particles(last, 0);
          particles(index, 1) = particles(last, 1);
          particles(index, 2) = particles(last, 2);
          particles(index, 3) = particles(last, 3);
          particles(index, 4) = particles(last, 4);
          particles(index, 5) = particles(last, 5);
        }

        --local_num;
      }
    }

    MPI_Allreduce(
      &local_num, &total_num, 1, MPI_INT, MPI_SUM, bunch.get_comm());
  }
}

void
populate_6d(Distribution& dist,
            Bunch& bunch,
            const_karray1d means,
            const_karray2d_row covariances)
{
  if (bunch.size() != bunch.get_local_num()) {
    throw std::runtime_error(
      "populate_6d: "
      "cannot populate bunches that has already lost particles.");
  }

  karray1d limits("limits", 6);
  for (int i = 0; i < 6; ++i) limits[i] = 0.0;

  populate_6d_truncated(dist, bunch, means, covariances, limits);
}

void
populate_6d_truncated(Distribution& dist,
                      Bunch& bunch,
                      const_karray1d means,
                      const_karray2d_row covariances,
                      const_karray1d limits)
{
#if 0
    multi_array_assert_size(means, 6, "populate_6d: means");
    multi_array_assert_size(covariances, 6, 6, "populate_6d: covariances");
    multi_array_assert_size(limits, 6, "populate_6d: limits");
#endif

  // deep copy from device to host
  bunch.checkout_particles();

  auto particles = bunch.get_host_particles();

  karray2d_row unit_covariances("unit_covariances", 6, 6);
  karray1d zero_means("zero_means", 6);

  bool truncated(false);

  for (int i = 0; i < 6; ++i) {
    double n = limits[i];

    if (n > 0) {
      truncated = true;

      double cutoff_integral =
        (exp(-n * n / 2.0)) *
        (sqrt(pi) * exp(n * n / 2.0) * erf(n / sqrt(2.0)) - sqrt(2.0) * n) /
        (sqrt(pi));

      unit_covariances(i, i) = 1.0 / (cutoff_integral * cutoff_integral);
    } else {
      unit_covariances(i, i) = 1.0;
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
          "populate_6d_truncated: "
          "maximum number of truncation iterations exceeded. "
          "Algorithm known to fail ~< 2.5 sigma.");
      }

      fill_unit_6d(dist, particles, unit_covariances, local_num, end);
      adjust_moments(bunch, zero_means, unit_covariances);
      strip_unit_6d(bunch, limits, total_num, local_num);
    }
  }

  adjust_moments(bunch, means, covariances);

  // copy to device
  bunch.checkin_particles();

  // check
  bunch.check_pz2_positive();
}

karray2d_row
get_correlation_matrix(const_karray2d_row one_turn_map,
                       double arms,
                       double brms,
                       double crms,
                       double beta,
                       std::array<int, 3> const& rms_index)
{

  for (int idx : rms_index) {
    if (idx < 0 || idx > 5)
      throw std::runtime_error("only valid rms indices (x=0, xp=1, y=2, "
                               "yp=3, z=4, dpp=5) are allowed");
  }

  karray2d_row correlation_matrix("corr_matrix", 6, 6);

  get_correlation_matrix_host(correlation_matrix.data(),
                              one_turn_map.data(),
                              arms,
                              brms,
                              crms,
                              beta,
                              rms_index);

  return correlation_matrix;
}

void
populate_transverse_gaussian(Distribution& dist,
                             Bunch& bunch,
                             const_karray1d means,
                             const_karray2d_row covariances,
                             double cdt)
{
  // deep copy from device to host
  bunch.checkout_particles();
  auto particles = bunch.get_host_particles();

  for (int i = 0; i < 4; ++i) {
    for (int p = 0; p < bunch.size(); ++p)
      particles(p, i) = dist.get_unit_gaussian();
  }

  // cdt
  for (int p = 0; p < bunch.size(); ++p)
    particles(p, 4) = dist.get_uniform(0.0, 1.0);

  // dpop
  for (int p = 0; p < bunch.size(); ++p)
    particles(p, 5) = dist.get_unit_gaussian();

  // copy of original means and covariances
  karray1d means_modified("means_modified");
  karray2d_row covariances_modified("covariances_modified");

  Kokkos::deep_copy(means_modified, means);
  Kokkos::deep_copy(covariances_modified, covariances);

  means_modified[4] = 0.0;

  // Symmetry requires no correlations with the cdt coordinate. Make a copy
  // of the covariance matrix and manually set all correlations to zero.
  for (int k = 0; k < 6; ++k)
    covariances_modified(k, 4) = covariances_modified(4, k) = 0.0;

  covariances_modified(4, 4) = cdt * cdt / 12.0;
  adjust_moments(bunch, means_modified, covariances_modified);

  // copy to device
  bunch.checkin_particles();
}

#if 0
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
populate_transverse_KV_GaussLong(Distribution &dist, Bunch &bunch, double epsilMax_x,
        double alpha_x, double beta_x, double epsilMax_y, double alpha_y, double beta_y,
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
// Pick first two points within a 2D circle, flat distribution
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
      // 4.0 is factor for x rms of uniform circular distribution
      const double a2X = std::sqrt((x1*x1 + x2*x2) * 4.0 * epsilMax_x*beta_x); // The amplitude.. X physical plane
      particles[part][Bunch::x] = a2X*std::sin(phi2X);
      particles[part][Bunch::xp] = (1.0/beta_x)*(a2X*std::cos(phi2X) - alpha_x*particles[part][Bunch::x]);
      //  Repeat in y,y' plane
      const double phi2Y = std::atan2(w,z);
      const double a2Y = std::sqrt((w*w + z*z) * 4.0 * epsilMax_y*beta_y); // The amplitude..
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



void
populate_longitudinal_boxcar(Distribution &dist, Bunch &bunch,   Const_MArray2d_ref one_turn_map, double length)
{
/// the corelation between the longitudinal and the transverse plane is neglected

    double cosmu=0.5*(one_turn_map[Bunch::cdt][Bunch::cdt]+one_turn_map[Bunch::dpop][Bunch::dpop]);
    if (fabs(cosmu)>1.) throw std::runtime_error("longitudinal modes: cosmu larger than zero");
    double sinmu=sqrt(1.-cosmu*cosmu);
    double alpha_z=(one_turn_map[Bunch::cdt][Bunch::cdt]-cosmu)/sinmu;
    double beta_z=one_turn_map[Bunch::cdt][Bunch::dpop]/sinmu;

    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];
    dist.fill_uniform(particles[boost::indices[range()][Bunch::dpop]], 0.0, 2.0*mconstants::pi);
    dist.fill_uniform(particles[boost::indices[range()][Bunch::cdt]], 0.0, 1.0);
    for (int part = 0; part < num_part; ++part) {
        double radius=length*sqrt(particles[part][Bunch::cdt]*(2.0-particles[part][Bunch::cdt]));
        double phase=particles[part][Bunch::dpop];
        particles[part][Bunch::cdt]=radius*cos(phase);
        particles[part][Bunch::dpop]=-radius*(sin(phase)+alpha_z*cos(phase))/beta_z;
    }

}

void
populate_longitudinal_uniform(Distribution &dist, Bunch &bunch,   double length)
{
    double half_length=0.5*length/bunch.get_reference_particle().get_beta();
    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];
    dist.fill_uniform(particles[boost::indices[range()][Bunch::cdt]], -half_length, half_length);

}

void
populate_transverseKV_logitudinalGaussian(Distribution &dist, Bunch &bunch,   Const_MArray2d_ref one_turn_map,
                             double radiusx,  double radiusy,    double ctrms)

{

// generates transversally a beam with radii radiusx and radiusy
// the transverse standard deviations will be xrms=radiusx/2,  yrms=radiusy/2
// valid for uncoupled linear maps

   double  cosmu=0.5*(one_turn_map[0][0]+one_turn_map[1][1]);
   if (fabs(cosmu)>1.) throw std::runtime_error("populate KV alpha_x: cosmu larger than zero");
   double sinmu=sqrt(1.-cosmu*cosmu);
   double alpha_x=(one_turn_map[0][0]-cosmu)/sinmu;
   double   beta_x=one_turn_map[0][1]/sinmu;

   cosmu=0.5*(one_turn_map[2][2]+one_turn_map[3][3]);
   if (fabs(cosmu)>1.) throw std::runtime_error("populate KV alpha_y: cosmu larger than zero");
   sinmu=sqrt(1.-cosmu*cosmu);
   double alpha_y=(one_turn_map[2][2]-cosmu)/sinmu;
   double beta_y=one_turn_map[2][3]/sinmu;

   cosmu=0.5*(one_turn_map[4][4]+one_turn_map[5][5]);
   if (fabs(cosmu)>1.) throw std::runtime_error("populate KV alpha_z: cosmu larger than zero");
   sinmu=sqrt(1.-cosmu*cosmu);
   double alpha_z=(one_turn_map[4][4]-cosmu)/sinmu;
   double beta_z=one_turn_map[4][5]/sinmu;



   MArray2d_ref particles(bunch.get_local_particles());
   int num_part = particles.shape()[0];
// generate 6 gaussian distributions
   dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::x]]);
   dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::xp]]);
   dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::y]]);
   dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::yp]]);
   dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::z]]);
   dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::zp]]);
   for (int part = 0; part < num_part; ++part) {
      double r=sqrt( particles[part][Bunch::x]*particles[part][Bunch::x]+
             particles[part][Bunch::xp]*particles[part][Bunch::xp]+
             particles[part][Bunch::y]*particles[part][Bunch::y]+
             particles[part][Bunch::yp]*particles[part][Bunch::yp]   );
      //transverse distribution
      if (r>0.){
        particles[part][Bunch::x] *= radiusx/r;
        particles[part][Bunch::xp] =
              (radiusx*particles[part][Bunch::xp]/r-alpha_x*particles[part][Bunch::x])/beta_x;


        particles[part][Bunch::y] *= radiusy/r;
        particles[part][Bunch::yp]=
          (radiusy*particles[part][Bunch::yp]/r-alpha_y*particles[part][Bunch::y])/beta_y;
      }
      else{ // very unlikely  case, what should be done???
        throw std::runtime_error(
          "generating KV fails, zero radius point for 4d sphere, need to be fixed, try different seed");
      }
     // longitudinal distribution

     particles[part][Bunch::z]*= ctrms;
     particles[part][Bunch::zp]=
           (ctrms*particles[part][Bunch::zp]-alpha_z*particles[part][Bunch::z])/beta_z;
  } //for part



}
#endif
