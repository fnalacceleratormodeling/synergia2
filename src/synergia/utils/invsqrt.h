#ifndef INVSQRT_H
#define INVSQRT_H

// invsqrt is an instrinsic in some math implementations
// we make it available here in preparation for future use

template<typename T >
    inline T invsqrt(T const& x) { return 1.0 / sqrt(x); }

#endif // INVSQRT_H

