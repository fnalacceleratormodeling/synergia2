#ifndef FLOATING_POINT_H_
#define FLOATING_POINT_H_

#include <cmath>

/// check a = b within the given tolerance
inline bool
floating_point_equal(double a, double b, double tolerance)
{
    if (std::abs(a) < tolerance) {
        return (std::abs(a - b) < tolerance);
    } else {
        return (std::abs((a - b) / a) < tolerance);
    }
}

/// check a <= b within the given tolerance
inline bool
floating_point_leq(double a, double b, double tolerance)
{
    return ((b - a) > -tolerance);
}

#endif /* FLOATING_POINT_H_ */
