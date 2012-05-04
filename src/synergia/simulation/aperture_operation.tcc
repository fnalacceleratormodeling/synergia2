#ifndef APERTURE_OPERATION_TCC_
#define APERTURE_OPERATION_TCC_

#include "synergia/foundation/math_constants.h"
#include <boost/math/special_functions/fpclassify.hpp>
#include <cmath>

template<typename T>
    void
    Aperture_operation::apply_impl(T & t, Bunch & bunch, int verbosity, Logger & logger)
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
//                        std::cout << "lost: " << part
//                                << "  " << particles[part][Bunch::x]
//                                << "  " << particles[part][Bunch::xp]
//                                << "  " << particles[part][Bunch::y]
//                                << "  " << particles[part][Bunch::yp]
//                                << std::endl;
                        // Move the last particle into this newly empty position
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
        deposit_charge(discarded * bunch.get_real_num() / bunch.get_total_num());
        bunch.set_local_num(local_num);
        double t1 = MPI_Wtime();
        if (verbosity > 5) {
            logger << "Aperture_operation: type = " << get_aperture_type()
                    << ", time = " << std::fixed << std::setprecision(3) << t1
                    - t0 << "s" << std::endl;
        }

    }

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
        deposit_charge(discarded * bunch.get_real_num() / bunch.get_total_num());
        bunch.set_local_num(local_num);
        double t1 = MPI_Wtime();
        if (verbosity > 5) {
            logger << "Aperture_operation: type = " << get_aperture_type()
                    << ", time = " << std::fixed << std::setprecision(3) << t1
                    - t0 << "s" << std::endl;
        }
    }

inline bool
Finite_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    bool keep = true;
    if (!boost::math::isfinite(particles[part][0])) keep = false;
    if (!boost::math::isfinite(particles[part][1])) keep = false;
    if (!boost::math::isfinite(particles[part][2])) keep = false;
    if (!boost::math::isfinite(particles[part][3])) keep = false;
    if (!boost::math::isfinite(particles[part][4])) keep = false;
    if (!boost::math::isfinite(particles[part][5])) keep = false;
    // negative pz^2 will give rise to non-finite numbers in fixed-t frames
    double p_scaled = 1.0 + particles[part][Bunch::dpop];
    double px_scaled = particles[part][Bunch::xp];
    double py_scaled = particles[part][Bunch::yp];
    double pz2_scaled = p_scaled * p_scaled - px_scaled * px_scaled - py_scaled
            * py_scaled;
    if (pz2_scaled < 0.0) {
        keep = false;
    }
    return !keep;
}

inline bool
Circular_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double r2 = particles[part][Bunch::x] * particles[part][Bunch::x]
            + particles[part][Bunch::y] * particles[part][Bunch::y];
    return (r2 > radius2);
}

inline bool
Elliptical_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    double scaled_r2 = particles[part][Bunch::x] * particles[part][Bunch::x]
            / h2 + particles[part][Bunch::y] * particles[part][Bunch::y] / v2;
    return (scaled_r2 > 1.0);
}

inline bool
Rectangular_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    return ((std::abs(particles[part][Bunch::x]) > 0.5 * width) || (std::abs(
            particles[part][Bunch::y]) > 0.5 * height));
}

inline bool
Polygon_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    std::complex<double > u(particles[part][Bunch::x],
            particles[part][Bunch::y]);
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
    return (theta_sum / (2.0 * mconstants::pi) < tiny);
}

inline bool
Wire_elliptical_aperture_operation::operator()(MArray2d_ref & particles,
        int part)
{
    double x = particles[part][Bunch::x];
    double y = particles[part][Bunch::y];
    double scaled_r2 = x * x / h2 + y * y / v2;
    bool retval;
    if (wire_x > 0.0) {
        retval = (scaled_r2 > 1.0) || ((x >= wire_x)
            && (x <= wire_x + wire_width)) || (x >= wire_x + wire_width + gap);
    } else if (wire_x < 0.0) {
        retval = (scaled_r2 > 1.0) || ((x <= wire_x)
            && (x >= wire_x - wire_width)) || (x <= wire_x - wire_width - gap);
    } else {
        throw std::runtime_error(
                "wire_elliptical_aperture_operation: wire_x and gap should not be zero");
    }
    return retval;
}

inline bool
Lambertson_aperture_operation::operator()(MArray2d_ref & particles, int part)
{
    bool retval;
    if (radius > 0.0) {
        retval = (particles[part][Bunch::x] >= radius);
    } else if (radius < 0.0) {
        retval = (particles[part][Bunch::x] <= radius);
    } else {
        throw std::runtime_error(
                "lambertson_aperture_operation: lambertson_aperture_radius should not be zero");
    }
    return retval;
}

#endif /* APERTURE_OPERATION_TCC_ */
