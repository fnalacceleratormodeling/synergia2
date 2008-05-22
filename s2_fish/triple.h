// -*- C++ -*-
#ifndef HAVE_TRIPLET_H
#define HAVE_TRIPLET_H

#include <vector>

template<class T>
class Triple
{
private:
    T storage[3];
public:
    Triple() { ; }

    Triple(T x, T y, T z) {
        storage[0] = x;
        storage[1] = y;
        storage[2] = z;
    }

    Triple(const Triple<T> &t) {
        storage[0] = t[0];
        storage[1] = t[1];
        storage[2] = t[2];
    }

    Triple(const T array[3]) {
        storage[0] = array[0];
        storage[1] = array[1];
        storage[2] = array[2];
    }

    Triple(const std::vector<T> &v) {
        storage[0] = v[0];
        storage[1] = v[1];
        storage[2] = v[2];
    }

    T operator[](int i) const {
        return storage[i];
    }

    T & operator[](int i) {
        return storage[i];
    }

    T get(int i) const {
        return storage[i];
    }

    T set(int i, T x) {
        return storage[i] = x;
    }

    operator T* () {
        return storage;
    }

    T* c_array() {
        return storage;
    }

    std::vector<T> vector() const {
        std::vector<T> retval(3);
        retval[0] = storage[0];
        retval[1] = storage[1];
        retval[2] = storage[2];
        return retval;
    }

    void scale(T constant) {
        storage[0] *= constant;
        storage[1] *= constant;
        storage[2] *= constant;
    }

    void add(T constant) {
        storage[0] += constant;
        storage[1] += constant;
        storage[2] += constant;
    }

    void add(Triple<T> triple) {
        storage[0] += triple[0];
        storage[1] += triple[1];
        storage[2] += triple[2];
    }
};

typedef Triple<int> Int3;
typedef Triple<double> Double3;
#endif
