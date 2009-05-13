#include "constraints.h"
#include <cmath>
#include <mpi.h>
#include "math_constants.h"

void
apply_longitudinal_periodicity_t(Macro_bunch_store &mbs)
{
    Array_1d<double> z = mbs.local_particles.slice(vector2(Range(4), Range()));
    for (Array_1d<double>::Iterator it = z.begin();
            it != z.end();
            ++it) {
        double tmp = *it + pi;
        if (tmp > 0) {
            *it = fmod(tmp, 2.0 * pi) - pi;
        } else {
            *it = fmod(tmp, 2.0 * pi) + pi;
        }
    }
}

void
apply_longitudinal_periodicity_z(Macro_bunch_store &mbs, double length)
{   
    double half_length=0.5*length;	
    Array_1d<double> z = mbs.local_particles.slice(vector2(Range(4), Range()));
    for (Array_1d<double>::Iterator it = z.begin();
            it != z.end();
            ++it) {
        double tmp = *it + half_length;
        if (tmp > 0) {
            *it = fmod(tmp, length) -half_length ;
        } else {
            *it = fmod(tmp, length) + half_length;
        }
    }
}



inline double sqr(double x)
{
    return x*x;
}

void
apply_circular_aperture(Macro_bunch_store &mbs, double radius)
{
    std::cout << "begin: local = " << mbs.local_num << ", total = " << mbs.total_num << std::endl;
    double scaled_radius2 = sqr(radius * mbs.units(0));
    //~ std::cout << "wtf: radius = " << radius << ", u(0) = " << mbs.units(0) << std::endl;
    int kept = 0;
    int discarded = 0;
    for (int i = 0; i < mbs.local_num; ++i) {
        bool try_discard = true;
        while (try_discard) {
            double r2 = sqr(mbs.local_particles(0, i)) +
                        sqr(mbs.local_particles(2, i));
            //~ std::cout << static_cast<int>(mbs.local_particles(6,i)) << ": " << scaled_radius2 << " ?< " << r2 << std::endl;
            if (r2 > scaled_radius2) {
                ++discarded;
                --mbs.local_num;
                if (i == mbs.local_num) {
                    // No more particles left
                    try_discard = false;
                } else {
                    // Move the last particle into this newly empty position
                    int last = mbs.local_num;
                    mbs.local_particles(0, i) = mbs.local_particles(0, last);
                    mbs.local_particles(1, i) = mbs.local_particles(1, last);
                    mbs.local_particles(2, i) = mbs.local_particles(2, last);
                    mbs.local_particles(3, i) = mbs.local_particles(3, last);
                    mbs.local_particles(4, i) = mbs.local_particles(4, last);
                    mbs.local_particles(5, i) = mbs.local_particles(5, last);
                    mbs.local_particles(6, i) = mbs.local_particles(6, last);
                }
            } else {
                ++kept;
                try_discard = false;
            }
        }
    }
    //~ std::cout << "kept = " << kept << ", discarded = " << discarded << std::endl;
    int rank, size;
    int old_total_num = mbs.total_num;
    rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Allreduce(&mbs.local_num, &mbs.total_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    } else {
        mbs.total_num = mbs.local_num;
    }
    std::cout << "proc " << rank << ": local = " << mbs.local_num
    << ", total = " << mbs.total_num << std::endl;
    mbs.total_current *= mbs.total_num / (1.0 * old_total_num);
}

