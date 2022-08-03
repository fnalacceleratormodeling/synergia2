#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/foundation/physical_constants.h"

int
main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  auto initargs =
    Kokkos::InitializationSettings{}; /* use default constructor */
  Kokkos::initialize(initargs);

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
      bunch_parts.access(0, 0) = 2.03125;
      bunch_parts.access(0, 2) = 2.03125;
      bunch_parts.access(0, 4) = 2.03125;
      // print intital coordinates
      logger << "before kick, particle at (x y z):" << '\n';
      for (int k = 0; k < 1; ++k) {
        logger << k << ": " << bunch_parts(k, 0) << ", " << bunch_parts(k, 2)
               << ", " << bunch_parts(k, 4) << '\n';
      }
      bunch.checkin_particles();
    }

    auto sc_ops = Space_charge_3d_open_hockney_options(64, 64, 64);
    sc_ops.comm_group_size = 1;
    sc_ops.green_fn = green_fn_t::linear;

    // space charge operator
    auto sc = Space_charge_3d_open_hockney(sc_ops);

    // set domain
    std::array<double, 3> offset = {2, 2, 2};
    std::array<double, 3> size = {4, 4, 4};
    sc.set_fixed_domain(offset, size);

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

      bunch.checkin_particles();
    }
  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
