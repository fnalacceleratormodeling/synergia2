#include "mpi.h"
#include <iostream>
#include <iomanip>
#ifdef USE_SIMPLE_TIMER_MEM
#include <fstream>
#include <string>
#endif // USE_SIMPLE_TIMER_MEM

inline double
simple_timer_current()
{
#ifdef USE_SIMPLE_TIMER
#ifdef USE_SIMPLE_TIMER_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#endif // USE_SIMPLE_TIMER_BARRIER
    return MPI_Wtime();
#else
    return 0.0;
#endif // USE_SIMPLE_TIMER
}

inline double
simple_timer_show(double t0, const char * label)
{
#ifdef USE_SIMPLE_TIMER
#ifdef USE_SIMPLE_TIMER_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#endif // USE_SIMPLE_TIMER_BARRIER
    double t1 = MPI_Wtime();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::ios_base::fmtflags old_flags(std::cout.flags());
        std::cout << "simple_timer:" << label << ":";
        std::cout << std::scientific << std::setprecision(8);
        std::cout << (t1-t0) << std::endl;
        std::cout.flags(old_flags);
    }
#ifdef USE_SIMPLE_TIMER_MEM
    // the order of fields in this file is vmsize, vmrss, "shared", "code",
    // "library", "data-stack", "dirty", units are 4k pages.
    // True for at least slf5,6 and ubuntu 10.04-11.10 
    std::ifstream statm("/proc/self/statm", std::ifstream::in);
    std::string vmsize, vmrss, vmstk, throwaway1, throwaway2,throwaway3;
    statm >> vmsize;
    statm >> vmrss;
    statm >> throwaway1; // shared
    statm >> throwaway2; // code;
    statm >> throwaway3; // library;
    statm >> vmstk; // data-stack
    statm.close();
    if (rank == 0) {
	std::cout << "simple_timer_mem:" << label << ":" << vmsize << ":" << vmrss << ":" << vmstk << std::endl;
    }
#endif // USE_SIMPLE_TIMER_MEM
    return MPI_Wtime();
#else
    return 0.0;
#endif // USE_SIMPLE_TIMER
}
