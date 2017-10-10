#ifndef APERTURE_OPERATION_TCC_
#define APERTURE_OPERATION_TCC_

#include "synergia/foundation/math_constants.h"
#include <boost/math/special_functions/fpclassify.hpp>
#include <cmath>
#include "synergia/utils/synergia_omp.h"

template<typename T>
    void
    Aperture_operation::apply_impl(T & t, Bunch & bunch, int verbosity, Logger & logger)
    {
        double t0 = MPI_Wtime();

        int nt;
        #pragma omp parallel
        { nt = omp_get_num_threads(); }

        MArray2d_ref particles(bunch.get_local_particles());
        int npart = bunch.get_local_num();

        int * discard = new int[npart];
        int * discard_count = new int[nt];

        int part_per_thread = npart / nt;
        #pragma omp parallel shared(nt, npart, particles, discard, discard_count)
        {
            int it = omp_get_thread_num();
            discard_count[it] = 0;

            int s = it * part_per_thread;
            int e = (it==nt-1) ? npart : (s+part_per_thread);

            for (int part = s; part < e; ++part)
            {
                if (t(particles, part)) 
                {
                    discard[part] = 1;
                    ++discard_count[it];
                }
                else
                {
                    discard[part] = 0;
                }
            }
        }

        int discarded = 0;
        for (int i=0; i<nt; ++i) discarded += discard_count[i];

        int head = 0;
        int tail = npart - 1;

        do
        {
            while (!discard[head] && head<tail) ++head;
            if (head >= tail) break;

            while ( discard[tail] && tail>head) --tail;
            if (head >= tail) break;

            particles[head][0] = particles[tail][0];
            particles[head][1] = particles[tail][1];
            particles[head][2] = particles[tail][2];
            particles[head][3] = particles[tail][3];
            particles[head][4] = particles[tail][4];
            particles[head][5] = particles[tail][5];
            particles[head][6] = particles[tail][6];

            ++head;
            --tail;

        } while(head < tail);

        double charge = (discarded > 0) ? discarded * bunch.get_real_num() / bunch.get_total_num() : 0.0;

        deposit_charge(charge);
        bunch.set_local_num(npart - discarded);

        double t1 = MPI_Wtime();
        if (verbosity > 5) 
        {
            logger << "Aperture_operation: type = " << get_aperture_type()
                   << ", discarded: " << discarded
                   << ", time = " << std::fixed << std::setprecision(3) << t1
                    - t0 << "s_n" << std::endl;
        }

        delete [] discard;
        delete [] discard_count;
    }


#if 0
template<typename T>
    void
    Aperture_operation::apply_impl(T & t, Bunch & bunch, int verbosity, Logger & logger)
    {
        double t0 = MPI_Wtime();

        int nt;
        #pragma omp parallel
        { nt = omp_get_num_threads(); }

        MArray2d_ref particles(bunch.get_local_particles());
        int npart = bunch.get_local_num();

        int * discard = new int[npart];
        int * discard_counts = new int[nt];

        int part_per_thread = npart / nt;
        #pragma omp parallel shared(nt, npart, particles, discard, discard_counts)
        {
            int it = omp_get_thread_num();

            int s = it * part_per_thread;
            int e = (it==nt-1) ? npart : (s+part_per_thread);

            discard_counts[it] = 0;

            for (int part = s; part < e; ++part)
            {
                if (t(particles, part)) 
                {
                    discard[s + discard_counts[it]] = part;
                    ++ discard_counts[it];
                }
            }
            //std::cout << "i = " << it << ", discarded = " << discard_counts[it] << "\n";

            #pragma omp barrier
        }

        int total_discarded = discard_counts[0];

        for (int t = 1; t < nt; ++t)
        {
            std::memcpy( discard + total_discarded, discard + t*part_per_thread,
	        discard_counts[t]*sizeof(int) );
            total_discarded += discard_counts[t];
        }

        for (int n = total_discarded - 1; n >= 0; --n)
        {
            // handle each particle in the list of discards

            if (discard[n] == npart-1) {
                // this is the last particle, just reduce count
                --npart;
            } else {
                // move the last particle into the position of this discarded particle then reduce count
                int idx = discard[n];
                int last = npart - 1;

                particles[idx][0] = particles[last][0];
                particles[idx][1] = particles[last][1];
                particles[idx][2] = particles[last][2];
                particles[idx][3] = particles[last][3];
                particles[idx][4] = particles[last][4];
                particles[idx][5] = particles[last][5];
                particles[idx][6] = particles[last][6];
                --npart;
            }
        }

        double charge = 0.0;
        if (total_discarded > 0) {
        	charge = total_discarded * bunch.get_real_num() / bunch.get_total_num();
        }
        deposit_charge(charge);
        bunch.set_local_num(npart);
        double t1 = MPI_Wtime();
        if (verbosity > 5) {
            logger << "Aperture_operation: type = " << get_aperture_type()
                   << ", discarded: " << total_discarded
                    << ", time = " << std::fixed << std::setprecision(3) << t1
                    - t0 << "s_n" << std::endl;
        }

        delete [] discard;
        delete [] discard_counts;
    }
#endif

template<typename T>
    void
    Aperture_operation::dump_particles(T & t, Bunch & bunch, int verbosity, Logger & logger)
    {
        double t0 = MPI_Wtime();
        MArray2d_ref particles(bunch.get_local_particles());
        int discarded = 0;
        int local_num = bunch.get_local_num();
        for (int part = 0; part < local_num; ++part) {
            bool try_discard = true;
            while (try_discard) {
                if (t(particles, part)) {
                    ++discarded;
                    --local_num;
                    if (part == local_num) {
                        // No more particles left
                        try_discard = false;
                    } else {
                        // Move the last particle into this newly empty position
                        std::cout << "extracted: " << part
                                << "  " << particles[part][Bunch::x]
                                << "  " << particles[part][Bunch::xp]
                                << "  " << particles[part][Bunch::y]
                                << "  " << particles[part][Bunch::yp]
                                << std::endl;
                        int last = local_num;
                        particles[part][0] = particles[last][0];
                        particles[part][1] = particles[last][1];
                        particles[part][2] = particles[last][2];
                        particles[part][3] = particles[last][3];
                        particles[part][4] = particles[last][4];
                        particles[part][5] = particles[last][5];
                        particles[part][6] = particles[last][6];
                    }
                } else {
                    try_discard = false;
                }
            }
        }
        double charge = 0.0;
        if (discarded > 0) {
        	charge = discarded * bunch.get_real_num() / bunch.get_total_num();
        }
        deposit_charge(charge);
        bunch.set_local_num(local_num);
        double t1 = MPI_Wtime();
        if (verbosity > 5) {
            logger << "Aperture_operation: type = " << get_aperture_type()
                    << ", time = " << std::fixed << std::setprecision(3) << t1
                    - t0 << "s_n" << std::endl;
        }
    }

inline bool
Finite_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    bool keep = true;
    double sum = particles[part][0] + particles[part][1] + particles[part][2]
            + particles[part][3] + particles[part][4] + particles[part][5];
    if (!boost::math::isfinite(sum)) {
        keep = false;
    }
    // negative pz^2 will give rise to non-finite numbers in fixed-t frames
    double p_scaled = 1.0 + particles[part][Bunch::dpop];
    double px_scaled = particles[part][Bunch::xp];
    double py_scaled = particles[part][Bunch::yp];
    double pz2_scaled = p_scaled * p_scaled - px_scaled * px_scaled
            - py_scaled * py_scaled;
    if (pz2_scaled < 0.0) {
        keep = false;
    }
    return !keep;
}

inline bool
Circular_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double xrel = particles[part][Bunch::x] - get_x_offset();
    double yrel = particles[part][Bunch::y] - get_y_offset();

    double r2 = xrel * xrel + yrel * yrel;
    return (r2 > radius2);
}

inline bool
Elliptical_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double xrel = particles[part][Bunch::x] - get_x_offset();
    double yrel = particles[part][Bunch::y] - get_y_offset();

    double scaled_r2 = xrel * xrel / h2 + yrel * yrel / v2;
    return (scaled_r2 > 1.0);
}

inline bool
Rectangular_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double xrel = particles[part][Bunch::x] - get_x_offset();
    double yrel = particles[part][Bunch::y] - get_y_offset();

    return ((std::abs(xrel) > 0.5 * width) || (std::abs(yrel) > 0.5 * height));
}

inline bool
Rectangular_with_ears_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double xrel = particles[part][Bunch::x] - get_x_offset();
    double yrel = particles[part][Bunch::y] - get_y_offset();

    // hopefully most particles will be within the rectangular area so check it first
    if ((std::abs(xrel) <= 0.5*width) &&
        (std::abs(yrel) <= 0.5*height))
        return 0;

    if ((std::abs(xrel) > 0.5*width+radius) ||
        (std::abs(yrel) > 0.5*height))
        return 1;

    if (std::abs(yrel) <= ear_offset)
        return (std::abs(xrel) > 0.5*width+radius);
 
    double xcirc = std::abs(xrel)-0.5*width;
    double ycirc = std::abs(yrel)-ear_offset;

    return ( xcirc*xcirc + ycirc*ycirc > radius*radius);
}

inline bool
Polygon_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double xrel = particles[part][Bunch::x] - get_x_offset();
    double yrel = particles[part][Bunch::y] - get_y_offset();
    double r2 = xrel * xrel + yrel * yrel;

    bool keep = true;
    if (r2 >= min_radius2) {
        std::complex<double > u(xrel, yrel);
        int index = 0;
        int size = vertices.size();
        double theta_sum = 0.0;
        while (index < size) {
            int index2 = index + 1;
            if (size == index2) index2 = 0;
            std::complex<double > v(vertices[index]);
            std::complex<double > w(vertices[index2]);
            double theta = arg((w - u) * conj(v - u));
            theta_sum += theta;
            ++index;
        }
        const double tiny = 1.0e-12;
        if (theta_sum / (2.0 * mconstants::pi) < tiny) keep = false;
    }
    return (!keep);
}

inline bool
Wire_elliptical_aperture_operation::operator()(MArray2d_ref & particles,
        int part)
{
    double xrel = particles[part][Bunch::x] - get_x_offset();
    double yrel = particles[part][Bunch::y] - get_y_offset();

    double scaled_r2 = xrel * xrel / h2 + yrel * yrel / v2;
    bool retval;
    if (wire_x > 0.0) {
        retval = (scaled_r2 > 1.0) || ((xrel >= wire_x)
            && (xrel <= wire_x + wire_width)) || (xrel >= wire_x + wire_width + gap);
    } else if (wire_x < 0.0) {
        retval = (scaled_r2 > 1.0) || ((xrel <= wire_x)
            && (xrel >= wire_x - wire_width)) || (xrel <= wire_x - wire_width - gap);
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_x and gap should not be zero");
    }
    return retval;
}

inline bool
Lambertson_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double xrel = particles[part][Bunch::x] - get_x_offset();

    bool retval;
    if (radius > 0.0) {
        retval = (xrel >= radius);
    } else if (radius < 0.0) {
        retval = (xrel <= radius);
    } else {
        throw std::runtime_error(
                "lambertson_aperture_operation: lambertson_aperture_radius should not be zero");
    }
    return retval;
}

#endif /* APERTURE_OPERATION_TCC_ */
