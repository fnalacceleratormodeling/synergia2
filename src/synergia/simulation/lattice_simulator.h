
#ifndef SYNERGIA_SIMULATION_LATTICE_SIMULATOR_H
#define SYNERGIA_SIMULATION_LATTICE_SIMULATOR_H

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice.h"
#include "synergia/libFF/ff_element.h"
#include "synergia/simulation/propagator.h"

#include "synergia/foundation/normal_form.h"
#include "synergia/foundation/trigon.h"

struct chromaticities_t {
  double momentum_compaction;

  double horizontal_chromaticity;
  double horizontal_chromaticity_prime;

  double vertical_chromaticity;
  double vertical_chromaticity_prime;

  double slip_factor;
  double slip_factor_prime;
};

namespace Lattice_simulator {
  constexpr const double default_closed_orbit_tolerance = 1.0e-13;

  void set_closed_orbit_tolerance(double tolerance);
  double get_closed_orbit_tolerance();

  // Both tune_linear_lattice() and tune_circular_lattice() set the frequency of
  // the rfcavities based on the momentum of the lattice reference particle.
  //
  // tune_linear_lattice uses the state of the reference particle as the
  // starting point for propagation.
  //
  // tune_circular_lattice() calculates a closed orbit and uses that as the
  // starting point.
  //
  // Both returns return an array of state, the transverse coordinates are the
  // final coordinates. cdt is the c * the total propagation time of the
  // particle which gives the total path length.

  // find the closed orbit for the circular lattice, propagate the lattice
  // reference particle through the lattice slices using the closed orbit, and
  // set the reference c*t for each lattice slice after tuning. return values is
  // the state for calcualted closed orbit note that all the rf cavities will be
  // set to 0 strength during the tuning process
  std::array<double, 6> tune_linear_lattice(Lattice& lattice);

  std::array<double, 6> tune_circular_lattice(Lattice& lattice);

  std::array<double, 6> tune_rfcavities(Lattice& lattice);

  // closed orbit
  std::array<double, 6> calculate_closed_orbit(Lattice const& lattice,
                                               double dpp = 0.0);

  // [tune_h, tune_v, c_delta_t]
  std::array<double, 3> calculate_tune_and_cdt(Lattice const& lattice,
                                               double dpp = 0.0);

  chromaticities_t get_chromaticities(Lattice const& lattice,
                                      double dpp = 1e-5);

  // get the full mapping of the one turn map
  template <unsigned int order = 2>
  TMapping<Trigon<double, order, 6>> get_one_turn_map(Lattice const& lattice,
                                                      double dpp = 0.0);

  // only the jacobian of the one turn map
  karray2d_row get_linear_one_turn_map(Lattice const& lattice);

  // [alpha, beta, psi]
  std::array<double, 3> map_to_twiss(karray2d_row map);

  template <unsigned int order>
  NormalForm<order> calculate_normal_form(Lattice const& lattice);

  double get_bucket_length(Lattice const& lattice);

  double get_rf_frequency(Lattice const& lattice);

  // the calculated lattice functions will be
  // written into the lattice elements
  //
  // overrides with Lattice& calculates the LF for elements
  // overrides with Propagator& do it for the element slices
  //
  void CourantSnyderLatticeFunctions(Lattice& lattice);

  void CourantSnyderLatticeFunctions(Propagator& prop);

  void calc_dispersions(Lattice& lattice);

  void calc_dispersions(Propagator& prop);

  // LF implementations
  template <class ELMS>
  void CourantSnyderLatticeFunctions_impl(Lattice& lattice, ELMS& elms);

  template <class ELMS>
  void calc_dispersions_impl(Lattice& lattice, ELMS& elms);

  // adjust tunes and chromaticities
  void adjust_tunes(Lattice& lattice,
                    double horizontal_tune,
                    double vertical_tune,
                    double tolerance = 1e-5);

  void adjust_chromaticities(Lattice& lattice,
                             double horizontal_chromaticity,
                             double vertical_chromaticity,
                             double tolerance = 1e-4,
                             int max_steps = 6);

}

#ifdef __CUDA_ARCH__

// no implementations for CUDA arch

#else

// implementations
namespace Lattice_simulator {
  template <unsigned int order>
  TMapping<Trigon<double, order, 6>>
  get_one_turn_map(Lattice const& lattice, double dpp)
  {
    using trigon_t = Trigon<double, order, 6>;

    // get the reference particle
    auto const& ref = lattice.get_reference_particle();

    // closed orbit
    auto probe = Lattice_simulator::calculate_closed_orbit(lattice, dpp);

    // comm world
    Commxx comm;

    // trigon bunch to get the one-turn-map
    bunch_t<trigon_t> tb(ref, comm.size(), comm);

    // design reference particle from the closed orbit
    auto ref_l = ref;
    ref_l.set_state(probe);
    tb.set_design_reference_particle(ref_l);

    auto tparts = tb.get_host_particles();

    // init value
    for (int i = 0; i < 6; ++i) tparts(0, i).set(probe[i], i);

    // check in
    tb.checkin_particles();

    // propagate trigon
    for (auto& ele : lattice.get_elements()) {
      if (ele.get_type() == element_type::rfcavity) {
        Lattice_element dup = ele;
        // dup.set_double_attribute("volt", 0.0);

        FF_element::apply(dup, tb);
      } else {
        FF_element::apply(ele, tb);
      }
    }

    // checkout particles
    tb.checkout_particles();

    // one-turn-map
    TMapping<trigon_t> map;
    for (int i = 0; i < trigon_t::dim; ++i) map[i] = tparts(0, i);

    return map;
  }

  template <unsigned int order>
  NormalForm<order>
  calculate_normal_form(Lattice const& lattice)
  {
    auto one_turn_map = get_one_turn_map<order>(lattice, 0.00);

    auto ref = lattice.get_reference_particle();
    double e0 = ref.get_total_energy();
    double pc0 = ref.get_momentum();
    double mass = ref.get_mass();

    NormalForm<order> nf(one_turn_map, e0, pc0, mass);
    return nf;
  }

  template <class ELMS>
  void
  CourantSnyderLatticeFunctions_impl(Lattice& lattice, ELMS& elms)
  {
    constexpr const int order = 2;
    using trigon_t = Trigon<double, order, 6>;

    const int ix = 0;
    const int ipx = 1;
    const int iy = 2;
    const int ipy = 3;

    auto const& ref = lattice.get_reference_particle();
    auto probe = calculate_closed_orbit(lattice);

    auto map = get_one_turn_map<order>(lattice);
    auto jac = map.jacobian();

    // .......... Check coupling ............................
    //::checkForCoupling(mtrx);

    // Calculate initial lattice functions ...
    // ... first horizontal

    double cs = (jac(ix, ix) + jac(ipx, ipx)) / 2.0;

    if (fabs(cs) > 1.0) {
      std::stringstream ss;
      ss << "*** ERROR ***                                     \n"
         << "*** ERROR *** LattSim::CourantSnyderLatticeFunctions \n"
         << "*** ERROR *** cos( psi_H ) = " << cs << "\n"
         << "*** ERROR *** Lattice is unstable.                \n"
         << "*** ERROR *** Cannot continue with calculation.   \n"
         << "*** ERROR ***                                     \n"
         << std::endl;

      throw std::runtime_error(ss.str());
    }

    double sn =
      (jac(ix, ipx) > 0.0) ? sqrt(1.0 - cs * cs) : -sqrt(1.0 - cs * cs);

    if (sn == 0.0) {
      std::stringstream ss;
      ss << "*** ERROR ***                                     \n"
         << "*** ERROR *** LattSim::CourantSnyderLatticeFunctions \n"
         << "*** ERROR *** Integer horizontal tune.            \n"
         << "*** ERROR ***                                     \n"
         << std::endl;

      throw std::runtime_error(ss.str());
    }

    double beta_x = jac(ix, ipx) / sn;
    double alpha_x = (jac(ix, ix) - jac(ipx, ipx)) / (2.0 * sn);

    // ... then vertical.
    cs = (jac(iy, iy) + jac(ipy, ipy)) / 2.0;

    if (fabs(cs) <= 1.0) {
      if (jac(iy, ipy) > 0.0)
        sn = sqrt(1.0 - cs * cs);
      else
        sn = -sqrt(1.0 - cs * cs);
    } else {
      std::stringstream ss;
      ss << "*** ERROR ***                                     \n"
         << "*** ERROR *** LattSim::CourantSnyderLatticeFunctions \n"
         << "*** ERROR *** cos( psi_V ) = " << cs << "\n"
         << "*** ERROR *** Lattice is unstable.                \n"
         << "*** ERROR *** Cannot continue with calculation.   \n"
         << "*** ERROR ***                                     \n"
         << std::endl;

      throw std::runtime_error(ss.str());
    }

    if (sn == 0.0) {
      std::stringstream ss;
      ss << "*** ERROR ***                                     \n"
         << "*** ERROR *** LattSim::CourantSnyderLatticeFunctions \n"
         << "*** ERROR *** Integer vertical tune.              \n"
         << "*** ERROR ***                                     \n"
         << std::endl;

      throw std::runtime_error(ss.str());
    }

    double beta_y = jac(iy, ipy) / sn;
    double alpha_y = (jac(iy, iy) - jac(ipy, ipy)) / (2.0 * sn);

    double beta0H = beta_x;
    double beta0V = beta_y;
    double alpha0H = alpha_x;
    double alpha0V = alpha_y;

    double oldpsiH = 0.0;
    double oldpsiV = 0.0;

    double tb = 0.0;
    double t = 0.0;
    double lng = 0.0;
    double psi_x = 0.0;
    double psi_y = 0.0;

    // trigon bunch
    Commxx comm;
    bunch_t<trigon_t> bunch(ref, comm.size(), comm);

    // design reference particle from the closed orbit
    auto ref_l = ref;
    ref_l.set_state(probe);
    bunch.set_design_reference_particle(ref_l);

    auto tparts = bunch.get_host_particles();

    // init value set to id map, ref points set to state
    // equivalent to jparticle.setState( particle.State() );
    for (int i = 0; i < 6; ++i) tparts(0, i).set(probe[i], i);

    // check in
    bunch.checkin_particles();

    // propagate trigon
    for (auto& elm : elms) {
      // At one time, dipoles with non-standard faces were discriminated
      // against and wouldn't have
      // their phase advance calculated.
      // bool is_regular =
      //     ( ( typeid(*lbe) != typeid(rbend)    ) &&
      //     (   typeid(*lbe) != typeid(CF_rbend) ) &&
      //     (   typeid(*lbe) != typeid(Slot)     ) &&
      //     (   typeid(*lbe) != typeid(srot)     ) &&
      //     (     (*lbe).hasStandardFaces()      )  );

      bool is_regular = true;

      // lng += elm.OrbitLength( particle );
      lng += elm.get_length();
      FF_element::apply(elm, bunch);

      auto mtrx = bunch.get_jacobian(0);

      tb = mtrx(ix, ix) * beta0H - mtrx(ix, ipx) * alpha0H;
      beta_x = (tb * tb + mtrx(ix, ipx) * mtrx(ix, ipx)) / beta0H;

      alpha_x = -1.0 *
                (tb * (mtrx(ipx, ix) * beta0H - mtrx(ipx, ipx) * alpha0H) +
                 mtrx(ix, ipx) * mtrx(ipx, ipx)) /
                beta0H;

      if (is_regular) {
        t = atan2(mtrx(ix, ipx), tb);

        // numerical round off errs introduce unphisical jumps in phase
        // while(t < oldpsiH) t += M_TWOPI;
        while (t < oldpsiH * (1. - 1.e-4)) t += mconstants::pi * 2;

        psi_x = oldpsiH = t;
      } else {
        psi_x = oldpsiH;
      }

      tb = mtrx(iy, iy) * beta0V - mtrx(iy, ipy) * alpha0V;
      beta_y = (tb * tb + mtrx(iy, ipy) * mtrx(iy, ipy)) / beta0V;

      alpha_y = -1.0 *
                (tb * (mtrx(ipy, iy) * beta0V - mtrx(ipy, ipy) * alpha0V) +
                 mtrx(iy, ipy) * mtrx(ipy, ipy)) /
                beta0V;

      if (is_regular) {
        t = atan2(mtrx(iy, ipy), tb);

        // numerical round off errs introduce unphisical jumps in phase
        // while(t < oldpsiV) t += M_TWOPI;
        while (t < oldpsiV * (1. - 1.e-4)) t += mconstants::pi * 2;

        psi_y = oldpsiV = t;
      } else {
        psi_y = oldpsiV;
      }

      elm.lf.arcLength = lng;
      elm.lf.beta.hor = beta_x;
      elm.lf.beta.ver = beta_y;
      elm.lf.alpha.hor = alpha_x;
      elm.lf.alpha.ver = alpha_y;
      elm.lf.psi.hor = psi_x;
      elm.lf.psi.ver = psi_y;

    } // End loop on lbe ...
  }

  template <class ELMS>
  void
  calc_dispersions_impl(Lattice& lattice, ELMS& elms)
  {
    const double dpp = 0.0005;

    constexpr const int order = 2;
    using trigon_t = Trigon<double, order, 6>;

    const int ix = 0;
    const int ipx = 1;
    const int iy = 2;
    const int ipy = 3;

    auto const& ref = lattice.get_reference_particle();

    // Preliminary steps ...
    auto probe1 = calculate_closed_orbit(lattice, 0.0);
    auto probe2 = calculate_closed_orbit(lattice, dpp);

    // propagate through elements
    Commxx comm;
    bunch_t<double> b1(ref, comm.size(), 1e9, comm);
    bunch_t<double> b2(ref, comm.size(), 1e9, comm);

    // design reference particle from the closed orbit
    auto ref1 = ref;
    ref1.set_state(probe1);
    b1.set_design_reference_particle(ref1);

    auto ref2 = ref;
    ref2.set_state(probe2);
    b2.set_design_reference_particle(ref2);

    // local particle data
    auto part1 = b1.get_host_particles();
    auto part2 = b2.get_host_particles();

    // init value
    for (int i = 0; i < 6; ++i) {
      part1(0, i) = probe1[i];
      part2(0, i) = probe2[i];
    }

    // check in
    b1.checkin_particles();
    b2.checkin_particles();

    // arcLength
    double lng = 0.0;

    // Attach initial dispersion data to the elements...
    for (auto& elm : elms) {
      lng += elm.get_length();

      FF_element::apply(elm, b1);
      FF_element::apply(elm, b2);

      b1.checkout_particles();
      b2.checkout_particles();

      std::array<double, 6> d;

      for (int i = 0; i < 6; ++i) d[i] = (part2(0, i) - part1(0, i)) / dpp;

      elm.lf.dispersion.hor = d[ix];
      elm.lf.dispersion.ver = d[iy];
      elm.lf.dPrime.hor = d[ipx];
      elm.lf.dPrime.ver = d[ipy];
      elm.lf.arcLength = lng;
    }
  }

}

#endif // __CUDA_ARCH

#endif // LATTICE_SIMULATOR_H
