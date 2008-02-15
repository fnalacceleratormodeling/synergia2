#ifndef HAVE_VECTOR_HELPER_H
#define HAVE_VECTOR_HELPER_H true

#include <vector>

template <typename T>
std::vector<T> vector1(T v0) {
    std::vector<T> retval(1);
    retval[0] = v0;
    return retval;
}

template <typename T>
std::vector<T> vector2(T v0, T v1) {
    std::vector<T> retval(2);
    retval[0] = v0;
    retval[1] = v1;
    return retval;
}

template <typename T>
std::vector<T> vector3(T v0, T v1, T v2) {
    std::vector<T> retval(3);
    retval[0] = v0;
    retval[1] = v1;
    retval[2] = v2;
    return retval;
}

template <typename T>
std::vector<T> vector4(T v0, T v1, T v2, T v3) {
    std::vector<T> retval(4);
    retval[0] = v0;
    retval[1] = v1;
    retval[2] = v2;
    retval[3] = v3;
    return retval;
}

#endif
