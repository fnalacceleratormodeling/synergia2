#include "mytimer.h"
#include <iostream>

#ifdef ENABLE_MYTIMER_CLOCK
#include <sys/times.h>
#endif

#ifdef ENABLE_MYTIMER_MPI
#include <mpi.h>
#endif

double timer_t0;
void
reset_timer()
{
    timer_t0 = time();
}

#ifdef ENABLE_MYTIMER_CLOCK
double
time()
{
    //~ timeb timestruct;
    //~ ftime(&timestruct);
    //~ return timestruct.time + timestruct.millitm / 1000.0;

    //~ return clock()/(1.0*CLOCKS_PER_SEC);

    tms my_tms;
    times(&my_tms);
    return my_tms.tms_stime / (1.0*sysconf(_SC_CLK_TCK));
}

void
timer(std::string message)
{
    double t1 = time();
    int rank;
    std::cout << "timer " << message << ": " << t1 - timer_t0 << std::endl;
}
#endif // ENABLE_MYTIMER_CLOCK

#ifdef ENABLE_MYTIMER_MPI
double
time()
{
    //~ timeb timestruct;
    //~ ftime(&timestruct);
    //~ return timestruct.time + timestruct.millitm / 1000.0;
    return MPI_Wtime();
}

void
timer(std::string message)
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
#endif // ENABLE_MYTIMER_MPI

#if !defined(ENABLE_MYTIMER_CLOCK) && !defined(ENABLE_MYTIMER_MPI)
double
time(){return 0.0; }

void
timer(std::string message) { return; }
#endif // neither defined
