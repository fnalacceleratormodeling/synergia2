#include "synergia/collective/space_charge_rectangular.h"
#include "synergia/foundation/physical_constants.h"

#include "synergia/utils/catch.hpp"

TEST_CASE("PointCharge", "[PointCharge]")
{

  auto logger = Logger(0, LoggerV::DEBUG);
  auto simlogger = Logger(0, LoggerV::INFO_STEP);

  const double mass = 100.0;
  const double total_energy = 125.0;
  const int total_num = 100;
  const double real_num = 2.0e12;

  /* scope to prevent deallocate after finalize errors for Kokkos views! */
  {
    Four_momentum fm(mass, total_energy);
    Reference_particle ref(pconstants::proton_charge, fm);

    auto bsim =
      Bunch_simulator::create_single_bunch_simulator(ref, 1, 1, Commxx());

    Bunch& bunch = bsim.get_bunch();
    {
      bunch.checkout_particles();
      auto bunch_parts = bunch.get_host_particles();
      bunch_parts.access(0, 0) = 0;
      bunch_parts.access(0, 2) = 0;
      bunch_parts.access(0, 4) = 0;
      // print intital coordinates
      logger << "before kick, particle at (x y z):" << '\n';
      for (int k = 0; k < 1; ++k) {
        logger << k << ": " << bunch_parts(k, 0) << ", " << bunch_parts(k, 2)
               << ", " << bunch_parts(k, 4) << '\n';
      }
      bunch.checkin_particles();
    }

    auto sc_ops = Space_charge_rectangular_options(
      std::array<int, 3>{17, 17, 17},
      std::array<double, 3>{4, 4, 4 * 0.59999999999999998});
    sc_ops.comm_group_size = 1;

    // space charge operator
    auto sc = Space_charge_rectangular(sc_ops);

    // apply space charge operator
    sc.apply(bsim, 1e-6, simlogger);

    {
      bunch.checkout_particles();
      auto bunch_parts = bunch.get_host_particles();
      // print final coordinates
      logger << "after kick, particle at (x y z):" << '\n';
      for (int k = 0; k < 1; ++k) {
        logger << k << ": " << bunch_parts(k, 0) << ", " << bunch_parts(k, 2)
               << ", " << bunch_parts(k, 4) << '\n';
      }
      logger << "after kick, particle momenta:" << '\n';
      for (int k = 0; k < 1; ++k) {
        logger << k << ": " << bunch_parts(k, 1) << ", " << bunch_parts(k, 3)
               << ", " << bunch_parts(k, 5) << '\n';
      }

      CHECK(bunch_parts(0, 1) == Approx(0).margin(.01));
      CHECK(bunch_parts(0, 3) == Approx(0).margin(.01));
      CHECK(bunch_parts(0, 5) == Approx(0).margin(.01));

      bunch.checkin_particles();
    }
  }
}
