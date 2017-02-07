// Generic SIMD Vector
#ifndef GSVECTOR_H_
#define GSVECTOR_H_

#include <boost/core/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>


#if 0
#undef GSV_SSE
#undef GSV_AVX
#undef GSV_V4D
#undef GSV_MIC
#endif


// helper
namespace detail
{
    template <class T, class E = void>
    struct VectorHelper
    { 
        static const size_t size() { return 1; }
        static T ld(const double *p) { return *p; } 
        static T st(double * p, const T & v) { *p = v; }
    };
}

// expression class
template <typename E, class T>
struct VecExpr
{
    typedef T vec_t;

    vec_t cal() const
    { return static_cast<E const&>(*this).cal(); } 

    operator E &()
    { return static_cast<      E &>(*this); }

    operator E const &() const
    { return reinterpret_cast<const E &>(*this); }
    //{ return static_cast<const E &>(*this); }
};

// the vector wrapper base class
template<class T>
struct Vec : public VecExpr<Vec<T>, T>
{
    T data;

    static const size_t size() { return detail::VectorHelper<T>::size(); }

    Vec(const double   d) : data( d ) { }
    Vec(const double * p) : data( detail::VectorHelper<T>::ld(p) ) { }

    void load (const double *p) { detail::VectorHelper<T>::ld(p); }
    void store(double *p) const { detail::VectorHelper<T>::st(p, data); }

    T & cal()       { return data; }
    T   cal() const { return data; }

    template <typename E>
    Vec(VecExpr<E, T> const & vec)
    {
        E const& v = vec;
        data = v.cal();
    }
};

// expression classes
template <typename E1, typename E2, class T, class E = void>
struct VecAdd : public VecExpr<VecAdd<E1, E2, T>, T>
{
    E1 const & _u; E2 const & _v;
    VecAdd(VecExpr<E1, T> const & u, VecExpr<E2, T> const & v) : _u(u), _v(v) { }
    typename VecExpr<VecAdd<E1, E2, T>, T>::vec_t cal() const { return _u.cal() + _v.cal(); }
};

template <typename E1, typename E2, class T>
struct VecSub : public VecExpr<VecSub<E1, E2, T>, T>
{
    E1 const & _u; E2 const & _v;
    VecSub(VecExpr<E1, T> const & u, VecExpr<E2, T> const & v) : _u(u), _v(v) { }
    typename VecExpr<VecSub<E1, E2, T>, T>::vec_t cal() const { return _u.cal() - _v.cal(); }
};

template <typename E1, typename E2, class T>
struct VecMul : public VecExpr<VecMul<E1, E2, T>, T>
{
    E1 const & _u; E2 const & _v;
    VecMul(VecExpr<E1, T> const & u, VecExpr<E2, T> const & v) : _u(u), _v(v) { }
    typename VecExpr<VecMul<E1, E2, T>, T>::vec_t cal() const { return _u.cal() * _v.cal(); }
};

template <typename E1, typename E2, class T>
struct VecDiv : public VecExpr<VecDiv<E1, E2, T>, T>
{
    E1 const & _u; E2 const & _v;
    VecDiv(VecExpr<E1, T> const & u, VecExpr<E2, T> const & v) : _u(u), _v(v) { }
    typename VecExpr<VecDiv<E1, E2, T>, T>::vec_t cal() const { return _u.cal() / _v.cal(); }
};

template <typename E1, typename E2, class T>
struct VecAddAssign : public VecExpr<VecAddAssign<E1, E2, T>, T>
{
    E1 & _u; E2 const & _v;
    VecAddAssign(VecExpr<E1, T> & u, VecExpr<E2, T> const & v) : _u(u), _v(v) { }
    typename VecExpr<VecAddAssign<E1, E2, T>, T>::vec_t cal() { _u.cal() = _u.cal() + _v.cal(); return _u; }
};

template <typename E1, typename E2, class T>
struct VecSubAssign : public VecExpr<VecSubAssign<E1, E2, T>, T>
{
    E1 & _u; E2 const & _v;
    VecSubAssign(VecExpr<E1, T> & u, VecExpr<E2, T> const & v) : _u(u), _v(v) { }
    typename VecExpr<VecSubAssign<E1, E2, T>, T>::vec_t cal() { _u = _u.cal() - _v.cal(); return _u; }
};

template <typename E, class T>
struct VecNeg : public VecExpr<VecNeg<E, T>, T>
{
    E const & _u;
    VecNeg(VecExpr<E, T> const & u) : _u(u) { }
    typename VecExpr<VecNeg<E, T>, T>::vec_t cal() const { return -(_u.cal()); }
};

template <typename E, class T>
struct VecSqrt : public VecExpr<VecSqrt<E, T>, T>
{
    E const & _u;
    VecSqrt(VecExpr<E, T> const & u) : _u(u) { }
    typename VecExpr<VecSqrt<E, T>, T>::vec_t cal() const { return sqrt(_u.cal()); }
};

// overload operators
//

template <typename E1, typename E2, class T>
inline
VecAdd<E1, E2, T> const
operator+ (VecExpr<E1, T> const & u, VecExpr<E2, T> const & v)
{ return VecAdd<E1, E2, T>(u, v); }

template <typename E1, typename E2, class T>
inline
VecSub<E1, E2, T> const
operator- (VecExpr<E1, T> const & u, VecExpr<E2, T> const & v)
{ return VecSub<E1, E2, T>(u, v); }

template <typename E1, typename E2, class T>
inline
VecMul<E1, E2, T> const
operator* (VecExpr<E1, T> const & u, VecExpr<E2, T> const & v)
{ return VecMul<E1, E2, T>(u, v); }

template <typename E1, typename E2, class T>
inline
VecDiv<E1, E2, T> const
operator/ (VecExpr<E1, T> const & u, VecExpr<E2, T> const & v)
{ return VecDiv<E1, E2, T>(u, v); }

template <typename E1, typename E2, class T>
inline
VecAddAssign<E1, E2, T> const
operator+= (VecExpr<E1, T>      & u, VecExpr<E2, T> const & v)
{ return VecAddAssign<E1, E2, T>(u, v); }

template <typename E1, typename E2, class T>
inline
VecSubAssign<E1, E2, T> const
operator-= (VecExpr<E1, T>      & u, VecExpr<E2, T> const & v)
{ return VecSubAssign<E1, E2, T>(u, v); }

template <typename E, class T>
inline
VecNeg<E, T> const
operator- (VecExpr<E, T> const & u)
{ return VecNeg<E, T>(u); }

template <typename E, class T>
inline
VecSqrt<E, T> const
sqrt (VecExpr<E, T> const & u)
{ return VecSqrt<E, T>(u); }



// specialization for different platforms
//

class Vec2d;
class Vec4d;
class vector4double;

//#include <x86intrin.h>
//#include <immintrin.h>

#if defined(GSV_SSE)
  #include "vectorclass.h"
#elif defined(GSV_AVX)
  #include "vectorclass.h"
#elif defined(GSV_QPX)
  #include <mass_simd.h>
#endif

namespace detail
{
    // specialization of helper class
    template <class T>
    struct VectorHelper<T, typename boost::enable_if<boost::is_same<T, Vec2d > >::type>
    { 
        static const size_t size() { return 2; }
        static T ld(const double *p) { T t; t.load_a(p); return t; }
        static T st(double * p, const T & v) { v.store_a(p); }
    };

    template <class T>
    struct VectorHelper<T, typename boost::enable_if<boost::is_same<T, Vec4d > >::type>
    { 
        static const size_t size() { return 4; }
        static T ld(const double *p) { T t; t.load_a(p); return t; }
        static T st(double * p, const T & v) { v.store_a(p); }
    };
}


// operations
template <typename E1, typename E2, class T>
struct VecAdd<E1, E2, T, typename boost::enable_if<boost::is_same<T, vector4double> >::type>
 : public VecExpr<VecAdd<E1, E2, T>, T>
{
    E1 const & _u; E2 const & _v;
    VecAdd(VecExpr<E1, T> const & u, VecExpr<E2, T> const & v) : _u(u), _v(v) { }
    typename VecExpr<VecAdd<E1, E2, T>, T>::vec_t cal() const { return vec_add(_u.cal(), _v.cal()); }
};

// stream operator
template <class T, typename boost::enable_if<boost::is_same<T, double> >::type>
inline std::ostream & operator << (std::ostream & out, Vec<T> & v)
{
    out << "(" << v.data << ")";
    return out;
}




// define the GSVector type

#if defined(GSV_SSE)
  typedef Vec<Vec2d >        GSVector;
#elif defined(GSV_AVX)
  typedef Vec<Vec4d >        GSVector;
#elif defined(GSV_QPX)
  typedef Vec<vector4double> GSVector;
#else
  typedef Vec<double>        GSVector;
#endif



#if 0

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

#endif

#endif  // GSVECTOR_H_
