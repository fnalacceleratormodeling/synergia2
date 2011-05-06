#include "mpi.h"
#include <iostream>

inline double
simple_timer_current()
{
#ifdef USE_SIMPLE_TIMER
  return MPI_Wtime();
#else
  return 0.0;
#endif // USE_SIMPLE_TIMER
}

inline double
simple_timer_show(double t0, const char * label)
{
#ifdef USE_SIMPLE_TIMER
  double t1 = MPI_Wtime();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << label << ":" << (t1-t0) << std::endl;
  }
  return MPI_Wtime();
#else
  return 0.0;
#endif // USE_SIMPLE_TIMER
}

