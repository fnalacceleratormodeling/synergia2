#ifndef UTILS_H_
#define UTILS_H_

#include <Kokkos_Core.hpp>
#include <iostream>
#include <mpi.h>
#include <vector>

#if defined BUILD_FD_SPACE_CHARGE_SOLVER
#include <petsc.h>
#include <string>
#endif

namespace synergia {

  /// Initialize synergia. This function initializes MPI, Kokkos and
  /// PETSc if it is included.
  /// @param argc number of arguments
  /// @param argv arguments given to the program
  static void
  initialize(int argc, char* argv[])
  {

#if defined BUILD_FD_SPACE_CHARGE_SOLVER
    auto ierr = PetscInitialize(
      &argc, &argv, (char*)0, std::string("synergia2-v3 program!\n").c_str());
    PetscCallAbort(PETSC_COMM_WORLD, PetscLogDefaultBegin());
#else
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      std::runtime_error("Could not initialize MPI!");
    }
#endif

    auto settings =
      Kokkos::InitializationSettings(); /* use default constructor */

    auto num_threads_chars = std::getenv("OMP_NUM_THREADS");
    if (num_threads_chars != nullptr) {
      settings.set_num_threads(std::stoi(num_threads_chars));
    }

    Kokkos::initialize(settings);

    return;
  }

  /// Finalize synergia. This function finalizes MPI, Kokkos and
  /// PETSc if it is included.
  static void
  finalize()
  {
    Kokkos::finalize();
#if defined BUILD_FD_SPACE_CHARGE_SOLVER
    PetscCallAbort(PETSC_COMM_WORLD, PetscLogView(PETSC_VIEWER_STDOUT_WORLD));
    auto ierr = PetscFinalize();
#else
    if (MPI_Finalize() != MPI_SUCCESS) {
      std::runtime_error("Could not finalize MPI!");
    }
#endif
    return;
  }

}

#endif
