// Generic SIMD Vector
#ifndef GSVECTOR_H_

#include "boost/multi_array.hpp"
typedef boost::multi_array_ref<double, 2> MArray2d_ref;

#if 0
#undef GSV_SSE
#undef GSV_AVX
#undef GSV_V4D
#undef GSV_MIC
#endif

#if defined(GSV_SSE)

#include "vectorclass.h"
class GSVector {
   public:
    static const int size = 2;
    static const int implementation = 1;
    Vec2d vec;
    GSVector(const MArray2d_ref& marray, size_t index0, size_t index1) : vec() {
        vec.load(&marray[index0][index1]);
    }
    GSVector(double val) : vec(val) {}
    GSVector(Vec2d const& vec) : vec(vec) {}
    GSVector(const double* array_start) : vec() { vec.load(array_start); }
    void store(MArray2d_ref& marray, size_t index0, size_t index1) const {
        vec.store(&marray[index0][index1]);
    }
    void store(double* array_start) const { vec.store(array_start); }
    void wtf() {
        std::cout << "GSVector::wtf " << vec[0] << "," << vec[1] << std::endl;
    }

    friend std::ostream & operator << (std::ostream & out, GSVector & v);
};
static inline GSVector operator+(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec + b.vec);
}
static inline GSVector operator+(GSVector const& a, double b) {
    return a + GSVector(b);
}
static inline GSVector operator+(double a, GSVector const& b) {
    return GSVector(a) + b;
}
static inline GSVector& operator+=(GSVector& a, GSVector const& b) {
    a.vec = a.vec + b.vec;
    return a;
}
static inline GSVector operator-(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec - b.vec);
}
static inline GSVector operator-(GSVector const& a, double b) {
    return GSVector(a.vec - b);
}
static inline GSVector operator-(double a, GSVector const& b) {
    return GSVector(a - b.vec);
}
static inline GSVector operator-(GSVector const& a) { return GSVector(-a.vec); }
static inline GSVector& operator-=(GSVector& a, GSVector const& b) {
    a.vec = a.vec - b.vec;
    return a;
}
static inline GSVector operator*(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec * b.vec);
}
static inline GSVector operator*(GSVector const& a, double b) {
    return GSVector(a.vec * b);
}
static inline GSVector operator*(double a, GSVector const& b) {
    return GSVector(a * b.vec);
}
static inline GSVector& operator*=(GSVector& a, GSVector const& b) {
    a.vec = a.vec * b.vec;
    return a;
}
static inline GSVector operator/(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec / b.vec);
}
static inline GSVector operator/(GSVector const& a, double b) {
    return GSVector(a.vec / b);
}
static inline GSVector operator/(double a, GSVector const& b) {
    return GSVector(a / b.vec);
}
static inline GSVector& operator/=(GSVector& a, GSVector const& b) {
    a.vec = a.vec / b.vec;
    return a;
}
static inline GSVector sqrt(GSVector const& a) { return GSVector(sqrt(a.vec)); }
static inline GSVector invsqrt(GSVector const& a) {
    return GSVector(1.0 / sqrt(a.vec));
}

inline std::ostream & operator << (std::ostream & out, GSVector & v)
{
    out << "(" << v.vec[0] << ", " << v.vec[1] << ")";
    return out;
}

#elif defined(GSV_AVX)

#include "vectorclass.h"
class GSVector {
public:
    static const int size = 4;
    static const int implementation = 2;
    Vec4d vec;
    GSVector(const MArray2d_ref& marray, size_t index0, size_t index1) : vec() {
        vec.load(&marray[index0][index1]);
    }
    GSVector(double val) : vec(val) {}
    GSVector(Vec4d const& vec) : vec(vec) {}
    GSVector(const double* array_start) : vec() { vec.load(array_start); }
    void store(MArray2d_ref& marray, size_t index0, size_t index1) const {
        vec.store(&marray[index0][index1]);
    }
    void store(double* array_start) const { vec.store(array_start); }
    void wtf() {
        std::cout << "GSVector::wtf " << vec[0] << "," << vec[1] << std::endl;
    }

    friend std::ostream & operator << (std::ostream & out, GSVector & v);
};
static inline GSVector operator+(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec + b.vec);
}
static inline GSVector operator+(GSVector const& a, double b) {
    return a + GSVector(b);
}
static inline GSVector operator+(double a, GSVector const& b) {
    return GSVector(a) + b;
}
static inline GSVector& operator+=(GSVector& a, GSVector const& b) {
    a.vec = a.vec + b.vec;
    return a;
}
static inline GSVector operator-(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec - b.vec);
}
static inline GSVector operator-(GSVector const& a, double b) {
    return GSVector(a.vec - b);
}
static inline GSVector operator-(double a, GSVector const& b) {
    return GSVector(a - b.vec);
}
static inline GSVector operator-(GSVector const& a) { return GSVector(-a.vec); }
static inline GSVector& operator-=(GSVector& a, GSVector const& b) {
    a.vec = a.vec - b.vec;
    return a;
}
static inline GSVector operator*(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec * b.vec);
}
static inline GSVector operator*(GSVector const& a, double b) {
    return GSVector(a.vec * b);
}
static inline GSVector operator*(double a, GSVector const& b) {
    return GSVector(a * b.vec);
}
static inline GSVector& operator*=(GSVector& a, GSVector const& b) {
    a.vec = a.vec * b.vec;
    return a;
}
static inline GSVector operator/(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec / b.vec);
}
static inline GSVector operator/(GSVector const& a, double b) {
    return GSVector(a.vec / b);
}
static inline GSVector operator/(double a, GSVector const& b) {
    return GSVector(a / b.vec);
}
static inline GSVector& operator/=(GSVector& a, GSVector const& b) {
    a.vec = a.vec / b.vec;
    return a;
}
static inline GSVector sqrt(GSVector const& a) { return GSVector(sqrt(a.vec)); }
static inline GSVector invsqrt(GSVector const& a) {
    return GSVector(1.0 / sqrt(a.vec));
}

inline std::ostream & operator << (std::ostream & out, GSVector & v)
{
    out << "(" << v.vec[0] << ", " << v.vec[1] << ", " << v.vec[2] << ", " << v.vec[3] << ")";
    return out;
}

#elif defined(GSV_V4D)

#include <mass_simd.h>
class GSVector {
public:
    static const int size = 4;
    static const int implementation = 1;
    vector4double vec;
//    GSVector(const MArray2d_ref& marray, size_t index0, size_t index1) : vec() {
//        vec.load(&marray[index0][index1]);
//    }
    GSVector(double val) : vec((vector4double)(val)) {}
    GSVector(vector4double const& vec) : vec(vec) {}
    // The following seems like it should be GSVector(const double * array_start),
    // but that does not compile.
    GSVector(double* array_start) : vec(vec_lda(0, array_start)) {}
//    void store(MArray2d_ref& marray, size_t index0, size_t index1) const {
//        vec.store(&marray[index0][index1]);
//    }
    void store(double* array_start) const { vec_sta(vec, 0, array_start); }
};
static inline GSVector operator+(GSVector const& a, GSVector const& b) {
    return GSVector(vec_add(a.vec, b.vec));
}
static inline GSVector operator+(GSVector const& a, double b) {
    return GSVector(vec_add(a.vec, (vector4double)(b)));
}
static inline GSVector operator+(double a, GSVector const& b) {
    return GSVector(vec_add((vector4double)(a), b.vec));
}
static inline GSVector& operator+=(GSVector& a, GSVector const& b) {
    a.vec = vec_add(a.vec, b.vec);
    return a;
}
static inline GSVector operator-(GSVector const& a, GSVector const& b) {
    return GSVector(vec_sub(a.vec, b.vec));
}
static inline GSVector operator-(GSVector const& a, double b) {
    return GSVector(vec_sub(a.vec, (vector4double)(b)));
}
static inline GSVector operator-(double a, GSVector const& b) {
    return GSVector(vec_sub((vector4double)(a), b.vec));
}
static inline GSVector operator-(GSVector const& a) { return GSVector(vec_neg(a.vec)); }
static inline GSVector& operator-=(GSVector& a, GSVector const& b) {
    a.vec = vec_sub(a.vec, b.vec);
    return a;
}
static inline GSVector operator*(GSVector const& a, GSVector const& b) {
    return GSVector(vec_mul(a.vec, b.vec));
}
static inline GSVector operator*(GSVector const& a, double b) {
    return GSVector(vec_mul(a.vec, (vector4double)(b)));
}
static inline GSVector operator*(double a, GSVector const& b) {
    return GSVector(vec_mul((vector4double)(a), b.vec));
}
static inline GSVector& operator*=(GSVector& a, GSVector const& b) {
    a.vec = vec_mul(a.vec, b.vec);
    return a;
}
static inline GSVector operator/(GSVector const& a, GSVector const& b) {
    return GSVector(vec_swdiv_nochk(a.vec, b.vec));
}
static inline GSVector operator/(GSVector const& a, double b) {
    return GSVector(vec_swdiv_nochk(a.vec, (vector4double)(b)));
}
static inline GSVector operator/(double a, GSVector const& b) {
    return GSVector(vec_swdiv_nochk((vector4double)(a), b.vec));
}
static inline GSVector& operator/=(GSVector& a, GSVector const& b) {
    a.vec = vec_swdiv_nochk(a.vec, b.vec);
    return a;
}
static inline GSVector sqrt(GSVector const& a) { return GSVector(sqrtd4(a.vec)); }
static inline GSVector invsqrt(GSVector const& a) { return GSVector(rsqrtd4(a.vec)); }

#elif defined(GSV_MIC)

#include <mic/micvec.h>
class GSVector {
public:
    static const int size = 8;
    static const int implementation = 4;
    F64vec8 vec;
    GSVector(MArray2d_ref& marray, size_t index0, size_t index1) : vec(&marray[index0][index1]) { }
    GSVector(double val) : vec(val) {}
    GSVector(F64vec8 const& vec) : vec(vec) {}
    GSVector(double* array_start) : vec(array_start) {  }
    void store(MArray2d_ref& marray, size_t index0, size_t index1) const 
    { _mm512_store_pd(&marray[index0][index1], vec); }
    void store(double* array_start) const 
    { _mm512_store_pd(array_start, vec); }
#if 0
    void wtf() {
        std::cout << "GSVector::wtf " << vec[0] << "," << vec[1] << std::endl;
    }
#endif

    friend std::ostream & operator << (std::ostream & out, GSVector & v);
};
static inline GSVector operator+(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec + b.vec);
}
static inline GSVector operator+(GSVector const& a, double b) {
    return a + GSVector(b);
}
static inline GSVector operator+(double a, GSVector const& b) {
    return GSVector(a) + b;
}
static inline GSVector& operator+=(GSVector& a, GSVector const& b) {
    a.vec = a.vec + b.vec;
    return a;
}
static inline GSVector operator-(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec - b.vec);
}
static inline GSVector operator-(GSVector const& a, double b) {
    return GSVector(a.vec - b);
}
static inline GSVector operator-(double a, GSVector const& b) {
    return GSVector(a - b.vec);
}
static inline GSVector operator-(GSVector const& a) { return GSVector(-a.vec); }
static inline GSVector& operator-=(GSVector& a, GSVector const& b) {
    a.vec = a.vec - b.vec;
    return a;
}
static inline GSVector operator*(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec * b.vec);
}
static inline GSVector operator*(GSVector const& a, double b) {
    return GSVector(a.vec * b);
}
static inline GSVector operator*(double a, GSVector const& b) {
    return GSVector(a * b.vec);
}
static inline GSVector& operator*=(GSVector& a, GSVector const& b) {
    a.vec = a.vec * b.vec;
    return a;
}
static inline GSVector operator/(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec / b.vec);
}
static inline GSVector operator/(GSVector const& a, double b) {
    return GSVector(a.vec / b);
}
static inline GSVector operator/(double a, GSVector const& b) {
    return GSVector(a / b.vec);
}
static inline GSVector& operator/=(GSVector& a, GSVector const& b) {
    a.vec = a.vec / b.vec;
    return a;
}
static inline GSVector sqrt(GSVector const& a) { return GSVector(sqrt(a.vec)); }
static inline GSVector invsqrt(GSVector const& a) {
    return GSVector(1.0 / sqrt(a.vec));
}

inline std::ostream & operator << (std::ostream & out, GSVector & v)
{
    //out << "(" << v.vec[0] << ", " << v.vec[1] << ", " << v.vec[2] << ", " << v.vec[3] << ")";
    return out;
}

#else  // no SIMD implementation defined

class GSVector {
   public:
    static const int size = 1;
    static const int implementation = 0;

    double vec;
    GSVector(const MArray2d_ref& marray, size_t index0, size_t index1)
        : vec(marray[index0][index1]) {}
    GSVector(double val) : vec(val) {}
    GSVector(const double* array_start) : vec(*array_start) {}
    void store(MArray2d_ref& marray, size_t index0, size_t index1) const {
        marray[index0][index1] = vec;
    }
    void store(double* array_start) const { *array_start = vec; }

    friend std::ostream & operator << (std::ostream & out, GSVector & v);
};
static inline GSVector operator+(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec + b.vec);
}
static inline GSVector operator+(GSVector const& a, double b) {
    return a + GSVector(b);
}
static inline GSVector operator+(double a, GSVector const& b) {
    return GSVector(a) + b;
}
static inline GSVector& operator+=(GSVector& a, GSVector const& b) {
    a.vec = a.vec + b.vec;
    return a;
}
static inline GSVector operator-(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec - b.vec);
}
static inline GSVector operator-(GSVector const& a, double b) {
    return GSVector(a.vec - b);
}
static inline GSVector operator-(double a, GSVector const& b) {
    return GSVector(a - b.vec);
}
static inline GSVector operator-(GSVector const& a) { return GSVector(-a.vec); }
static inline GSVector& operator-=(GSVector& a, GSVector const& b) {
    a.vec = a.vec - b.vec;
    return a;
}
static inline GSVector operator*(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec * b.vec);
}
static inline GSVector operator*(GSVector const& a, double b) {
    return GSVector(a.vec * b);
}
static inline GSVector operator*(double a, GSVector const& b) {
    return GSVector(a * b.vec);
}
static inline GSVector& operator*=(GSVector& a, GSVector const& b) {
    a.vec = a.vec * b.vec;
    return a;
}
static inline GSVector operator/(GSVector const& a, GSVector const& b) {
    return GSVector(a.vec / b.vec);
}
static inline GSVector operator/(GSVector const& a, double b) {
    return GSVector(a.vec / b);
}
static inline GSVector operator/(double a, GSVector const& b) {
    return GSVector(a / b.vec);
}
static inline GSVector& operator/=(GSVector& a, GSVector const& b) {
    a.vec = a.vec / b.vec;
    return a;
}
static inline GSVector sqrt(GSVector const& a) { return GSVector(sqrt(a.vec)); }
static inline GSVector invsqrt(GSVector const& a) {
    return GSVector(1.0 / sqrt(a.vec));
}

inline std::ostream & operator << (std::ostream & out, GSVector & v)
{
    out << "(" << v.vec << ")";
    return out;
}

#endif

#define GSVECTOR_H_
#endif  // GSVECTOR_H_
