#include "synergia/utils/catch.hpp"

#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/collective/tests/rod_bunch.h"

TEST_CASE("real_apply_full_lowgamma", "[Rod_bunch]")
{
  auto logger = Logger(0, LoggerV::DEBUG);
  auto simlogger = Logger(0, LoggerV::INFO_STEP);

  const int gridx = 64;
  const int gridy = 64;
  const int gridz = 64;

  const double time_fraction = 1.0;
  const double step_length = 0.1;

  // bunch
  Rod_bunch_fixture_lowgamma fixture;

  auto& bunch = fixture.bsim.get_bunch();
  auto const& ref = bunch.get_reference_particle();
  auto parts = bunch.get_host_particles();

  const double beta = ref.get_beta();
  const double gamma = ref.get_gamma();
  const double betagamma = beta * gamma;

  const double time_step = step_length / (beta * pconstants::c);
  // const double bunchlen = bunch.get_longitudinal_boundary().second;
  const double bunchlen = 0.1;

  fixture.bsim.get_bunch().write_file(std::string("particles_before_kick.h5"));

  // check out particles before print
  bunch.checkout_particles();

  // print intital coordinates
  logger << "real_apply_full_lowgamma first four particles (x y z):" << '\n';
  for (int k = 0; k < 4; ++k) {
    logger << k << ": " << parts(k, 0) << ", " << parts(k, 2) << ", "
           << parts(k, 4) << '\n';
  }

  logger << "last four particles (x y z):" << '\n';
  for (int k = bunch.get_local_num() - 4; k < bunch.get_local_num(); ++k) {
    logger << k << ": " << parts(k, 0) << ", " << parts(k, 2) << ", "
           << parts(k, 4) << '\n';
  }

  logger << '\n';

  // space charge operator
  auto sc_ops = Space_charge_3d_open_hockney_options(gridx, gridy, gridz);
  sc_ops.comm_group_size = 1;
  sc_ops.green_fn = green_fn_t::linear;

  auto sc = Space_charge_3d_open_hockney(sc_ops);

  // set domain
  std::array<double, 3> offset = {0, 0, 0};
  std::array<double, 3> size = {
    parts(0, 0) * 4, parts(0, 0) * 4, bunchlen / beta};
  sc.set_fixed_domain(offset, size);

  // apply space charge operator
  sc.apply(fixture.bsim, time_step, simlogger);

  // check out particles
  bunch.checkout_particles();

  // print
  logger << "after sc::apply : bunch.local_particles(0, 0): " << parts(0, 0)
         << '\n';

  // Rod of charge Q over length L
  // E field at radius r $$ E = \frac{1}{2 \pi \epsilon_0} \frac{Q}{L}
  // \frac{1}{r} $$ B field at radius r $$ E = \frac{\mu_0}{2 \pi } \frac{Q
  // v}{L} \frac{1}{r} $$ Net EM force on electric+magnetic on probe of charge q
  // from E-B cancellation
  // $$ F = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L} \frac{1}{r}
  // travel over distance D at velocity v
  // \frac{\Delta p}{p} = \frac{1}{2 \pi \epsilon_0 \gamma^2} \frac{qQ}{L}
  // \frac{D}{m v^2} \frac{1}{r} convert to usual units \frac{\Delta p}{p} =
  // \frac{2 N r_p}{L \beta^2 \gamma^3} \frac{D}{r}

  double L = bunchlen;
  double N = bunch.get_real_num();

  logger << "L: " << L << '\n';
  logger << "N: " << N << '\n';
  logger << "step_length: " << step_length << '\n';
  logger << "beta: " << beta << '\n';
  logger << "gamma: " << gamma << '\n';
  logger << "betagamma: " << betagamma << '\n';
  logger << "x: " << parts(0, Bunch::x) << '\n';

  double computed_dpop =
    ((2.0 * N * pconstants::rp) / (L * betagamma * betagamma * gamma)) *
    (step_length / parts(0, Bunch::x));

  logger << "computed dpop: " << computed_dpop << '\n';
  logger << "particle dpop: " << parts(0, 1) << '\n';

  CHECK(parts(0, Bunch::xp) == Approx(computed_dpop).margin(.01));

  int nkicks = 0;
  for (int k = 0; k < bunch.get_local_num(); ++k) {
    if ((parts(k, 1) != 0.0) || (parts(k, 3) != 0.0)) {
      ++nkicks;

      if (nkicks < 10) {
        logger << "kick: " << nkicks << ", particle " << k << ": "
               << parts(k, 0) << ", " << parts(k, 1) << ", " << parts(k, 2)
               << ", " << parts(k, 3) << ", " << parts(k, 4) << ", "
               << parts(k, 5) << '\n';
      }
    }
  }

  fixture.bsim.get_bunch().write_file(std::string("particles_after_kick.h5"));
}
