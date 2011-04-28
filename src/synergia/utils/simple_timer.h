#include "mpi.h"
#include <iostream>

inline void
simple_timer_reset(double & t)
{
  t = MPI_Wtime();
}

inline void
simple_timer_show(double & t, const char * label)
{
  double t1 = MPI_Wtime();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << label << ":" << (t1-t) << std::endl;
  }
  t = MPI_Wtime();
}
    
