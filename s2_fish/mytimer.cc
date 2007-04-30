#include "mytimer.h"
#include <iostream>
#include <mpi.h>

double
time()
{
    timeb timestruct;
    ftime(&timestruct);
    return timestruct.time + timestruct.millitm / 1000.0;
}

double timer_t0;
void
reset_timer()
{
    timer_t0 = time();
}

void
nottimer(std::string message)
{
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = time();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "timer " << message << ": " << t1 - timer_t0 << std::endl;
        timer_t0 = time();
    }
}

void
timer(std::string message) { return; }
