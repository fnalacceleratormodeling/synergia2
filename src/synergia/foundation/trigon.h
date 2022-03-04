#ifndef TRIGON_H
#define TRIGON_H

#include "synergia/foundation/trigon_traits.h"

#include <algorithm>
#include <complex>
#include <iostream> // jfa remove me!!!!
#include <unordered_map>

#include <Kokkos_Core.hpp>
//#include <Kokkos_UnorderedMap.hpp>

#include <Eigen/Eigen>

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/kokkos_types.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/json.h"

#include <cereal/cereal.hpp>
#include <cereal/types/complex.hpp>

template<class T, size_t SIZE>
struct arr_t
{
    T data_[SIZE] = {};

    KOKKOS_INLINE_FUNCTION
    constexpr size_t size() const
    { return SIZE; }

    KOKKOS_INLINE_FUNCTION
    void fill(T t)
    { for(size_t i=0; i<SIZE; ++i) data_[i] = t; }

    KOKKOS_INLINE_FUNCTION
    T& at(size_t idx) { return data_[idx]; }

    KOKKOS_INLINE_FUNCTION
    T const& at(size_t idx) const { return data_[idx]; }

    KOKKOS_INLINE_FUNCTION
    T& operator[](size_t idx) { return data_[idx]; }

    KOKKOS_INLINE_FUNCTION
    T const& operator[](size_t idx) const { return data_[idx]; }

    KOKKOS_INLINE_FUNCTION
    T*       begin()       { return data_; }

    KOKKOS_INLINE_FUNCTION
    T const* begin() const { return data_; }

    KOKKOS_INLINE_FUNCTION
    T*       end()         { return data_ + SIZE; }

    KOKKOS_INLINE_FUNCTION
    T const* end() const   { return data_ + SIZE; }

    // conversion
    template<class U>
    KOKKOS_INLINE_FUNCTION
    void from(arr_t<U, SIZE> const& o)
    { for(int i=0; i<SIZE; ++i) data_[i] = static_cast<T>(o.data_[i]); }

    template<class U>
    KOKKOS_INLINE_FUNCTION
    arr_t<U, SIZE> to() const
    { 
        arr_t<U, SIZE> ret;
        for(int i=0; i<SIZE; ++i) ret.data_[i] = static_cast<U>(data_[i]); 
        return ret;
    }

    template<class AR>
    void serialize(AR& ar)
    {
        ar(data_);
    }
};

template <typename T, unsigned int Power, unsigned int Dim>
class Trigon;

template <typename TRIGON>
class TMapping;

namespace trigon_impl
{
    template<class T, size_t SIZE>
    KOKKOS_INLINE_FUNCTION
    void sort(arr_t<T, SIZE> & arr)
    {
        for(size_t i=0; i<SIZE-1; ++i)
        {
            for(size_t j=i+1; j<SIZE; ++j)
            {
                if (arr[i] > arr[j])
                {
                    T temp = arr[i];
                    arr[i] = arr[j];
                    arr[j] = temp;
                }
            }
        }
    }

    //__constant__ double map[100];
}

template<class T, size_t N>
KOKKOS_INLINE_FUNCTION
bool operator==(arr_t<T, N> const& lhs, arr_t<T, N> const& rhs)
{
    for(size_t i=0; i<N; ++i)
        if (lhs[i] != rhs[i]) return false;
    return true;
}

KOKKOS_INLINE_FUNCTION
constexpr int
factorial(int n)
{
    return n <= 1 ? 1 : (n * factorial(n - 1));
}

template <typename T, unsigned int Power, unsigned int Dim>
class Trigon;

KOKKOS_INLINE_FUNCTION
constexpr unsigned int
array_length(unsigned int i)
{
    return (i == 0) ? 1 : i;
}

// exp array with the dim of power, where the elements are the index of the
// coordinates. E.g.,
// [1, 2 ,3] is the exp array for a 3rd order term for xyz
// [1, 1, 0] => x*x => x^2
// [1, 2, 2] => x*y*y => xy^2
template <unsigned int Power>
using Index_t = arr_t<size_t, array_length(Power)>;

// array holding the exp indices for each term in the terms array
template <unsigned int Power, unsigned int Dim>
using Indices_t = arr_t<Index_t<Power>, Trigon<double, Power, Dim>::count>;

template <unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<Power==0, Indices_t<Power, Dim>>
indices()
{
    Indices_t<Power, Dim> retval;
    retval[0][0] = 0;
    return retval;
}

template <unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<Power==1, Indices_t<Power, Dim>>
indices()
{
    Indices_t<Power, Dim> retval;
    for (size_t i=0; i<Dim; ++i) retval[i][0] = i;
    return retval;
}

#if 0
template <unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<((Power == 1) || (Power == 0)), Indices_t<Power, Dim>>
indices()
{
    Indices_t<Power, Dim> retval;
    if (Power == 0) {
        retval[0][0] = 0;
    } else {
        for (size_t i = 0; i < Dim; ++i) {
            retval[i][0] = i;
        }
    }

    return retval;
}
#endif

template <unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(Power > 1), Indices_t<Power, Dim>>
indices()
{
    Indices_t<Power - 1, Dim> subindices = indices<Power - 1, Dim>();
    Indices_t<Power, Dim> retval;
    size_t count = 0;
    for (size_t i = 0; i < Dim; ++i) {
        for (size_t j = 0; j < subindices.size(); ++j) {
            if (subindices[j][0] >= i) {
                arr_t<size_t, Power> entry;
                entry[0] = i;
                for (size_t k = 0; k < (Power - 1); ++k) {
                    entry[k + 1] = subindices[j][k];
                }
                retval[count] = entry;
                ++count;
            }
        }
    }

    return retval;
}

template <unsigned int Length>
struct Array_hash
{
    constexpr static size_t max_dim = 8;

    KOKKOS_INLINE_FUNCTION
    size_t operator()(arr_t<size_t, Length> const& arr) const
    {
        size_t sum = 0;
        size_t mult = 1;
        for (auto&& i : arr) {
            sum += mult * i;
            mult *= max_dim;
        }
        return sum;
    }
};

template <unsigned int Power, unsigned int Dim>
using Map_t = std::unordered_map<arr_t<size_t, Power>, size_t, 
        Array_hash<Power>>;

#if 0
template <unsigned int Power, unsigned int Dim>
using Map_t = Kokkos::UnorderedMap<arr_t<size_t, Power>, size_t, 
        Kokkos::DefaultExecutionSpace,
        Array_hash<Power>>;
#endif

template <unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
const Indices_t<Power, Dim>&
canonical_to_index()
{
    static const Indices_t<Power, Dim> ind = indices<Power, Dim>();
    return ind;
}

template <unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
const Map_t<Power, Dim>
fill_index_to_canonical()
{
    Map_t<Power, Dim> map;
    auto the_indices = canonical_to_index<Power, Dim>();
    //auto the_indices = indices<Power, Dim>();

    for (size_t i = 0; i < the_indices.size(); ++i) {
        map[the_indices[i]] = i;
    }

    return map;
}

template <unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
const Map_t<Power, Dim>&
index_to_canonical()
{
    static const Map_t<Power, Dim> map = fill_index_to_canonical<Power, Dim>();
    return map;
}

KOKKOS_INLINE_FUNCTION
double term_to_json_val(double const& term)
{ return term; }

KOKKOS_INLINE_FUNCTION
std::string term_to_json_val(std::complex<double> const& term)
{ return "(" + std::to_string(term.real()) + ',' + std::to_string(term.imag()) + ")"; }

template <typename T, unsigned int Power, unsigned int Dim>
class Trigon
{
public:

    using data_type = T;
    static constexpr unsigned int dim = Dim;

    static constexpr unsigned int count =
        factorial(Dim+Power-1) / (factorial(Dim-1) * factorial(Power));

    typedef arr_t<T, count> Terms_t;

    Trigon<T, Power - 1, Dim> lower;
    Terms_t terms;

    KOKKOS_INLINE_FUNCTION
    Trigon() : lower() 
    { 
        terms.fill(0); 
    }

    KOKKOS_INLINE_FUNCTION
    Trigon(T val) : lower(val) 
    { 
        terms.fill(0); 
    }

    KOKKOS_INLINE_FUNCTION
    Trigon(T val, size_t index) : lower(val)
    {
        terms.fill(0);
        get_subpower<1>().terms[index] = 1; // jfa fixme
    }

    // implicit conversion
    template<class U>
    KOKKOS_INLINE_FUNCTION
    explicit Trigon(Trigon<U, Power, Dim> const& o)
    : terms(), lower(o.lower) { terms.from(o.terms); }

    KOKKOS_INLINE_FUNCTION
    static constexpr unsigned int power() { return Power; }

    KOKKOS_INLINE_FUNCTION
    const T& value() const { return get_subpower<0>().terms[0]; }

    KOKKOS_INLINE_FUNCTION
    T& value() { return get_subpower<0>().terms[0]; }

    KOKKOS_INLINE_FUNCTION
    void set(T val)
    {
        terms.fill(0);
        lower.set(val);
    }

    KOKKOS_INLINE_FUNCTION
    void set(T val, size_t index)
    {
        set(val);
        get_subpower<1>().terms[index] = 1;
    }

    template <unsigned int Subpower>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if_t<(Subpower < Power), Trigon<T, Subpower, Dim>&>
    get_subpower()
    { return lower.template get_subpower<Subpower>(); }

    template <unsigned int Subpower>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if_t<(Subpower < Power), Trigon<T, Subpower, Dim> const&>
    get_subpower() const
    { return lower.template get_subpower<Subpower>(); }

    template <unsigned int Subpower>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if_t<(Subpower == Power), Trigon<T, Subpower, Dim>&>
    get_subpower()
    { return *this; }

    template <unsigned int Subpower>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if_t<(Subpower == Power), Trigon<T, Subpower, Dim> const&>
    get_subpower() const
    { return *this; }

    template <typename F>
    KOKKOS_INLINE_FUNCTION
    void each_term(F f)
    {
        lower.each_term(f);
        auto inds = indices<Power, Dim>();
        for(int i=0; i<terms.size(); ++i) f(i, inds[i], terms[i]);
    }

    KOKKOS_INLINE_FUNCTION
    void set_term(size_t idx, T const& val)
    { terms[idx] = val; }

    KOKKOS_INLINE_FUNCTION
    T get_term(size_t idx)
    { return terms[idx]; }

    KOKKOS_INLINE_FUNCTION
    void set_term(unsigned int power, size_t idx, T const& val)
    {
        if (power == Power) set_term(idx, val);
        else if (power == 0) get_subpower<0>().set_term(idx, val);
        else if (power < Power) lower.set_term(power, idx, val);
    }

    KOKKOS_INLINE_FUNCTION
    T get_term(unsigned int power, size_t idx)
    {
        if (power == Power) return get_term(idx);
        else if (power == 0) return get_subpower<0>().get_term(idx);
        else if (power < Power) return lower.get_term(power, idx);
        return T{};
    }

    // keep terms with power in [lower, upper]
    KOKKOS_INLINE_FUNCTION
    void filter(unsigned int p_lower, unsigned int p_upper)
    {
        if (Power>p_upper || Power<p_lower) terms.fill(0);
        lower.filter(p_lower, p_upper);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator== (T rhs) const
    { 
        double const eps = 1e5 * std::numeric_limits<double>::epsilon();
        if (abs(value() - rhs) > eps) return false;

        for(auto const& v : terms)
            if (abs(v-rhs) > eps) return false;

        return (lower == rhs);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator!= (T rhs) const
    { return !((*this) == rhs); }

    KOKKOS_INLINE_FUNCTION
    bool operator< (double rhs) const
    { return value() < rhs; }

    KOKKOS_INLINE_FUNCTION
    bool operator> (double rhs) const
    { return value() > rhs; }

    KOKKOS_INLINE_FUNCTION
    bool operator<= (double rhs) const
    { return value() <= rhs; }

    KOKKOS_INLINE_FUNCTION
    bool operator>= (double rhs) const
    { return value() >= rhs; }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim>& operator+=(Trigon<T, Power, Dim> const& t)
    {
        lower += t.lower;
        for (size_t i = 0; i < terms.size(); ++i) {
            terms[i] += t.terms[i];
        }
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim>& operator+=(T val)
    {
        lower += val;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator+(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval += t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator+(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval += val;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator-() const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval.lower = -lower;
        for (size_t i = 0; i < terms.size(); ++i) {
            retval.terms[i] = -terms[i];
        }
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim>& operator-=(Trigon<T, Power, Dim> const& t)
    {
        lower -= t.lower;
        for (size_t i = 0; i < terms.size(); ++i) {
            terms[i] -= t.terms[i];
        }
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim>& operator-=(T val)
    {
        lower -= val;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator-(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval -= t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator-(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval -= val;
        return retval;
    }

    template <unsigned int P1, unsigned int P2>
    KOKKOS_INLINE_FUNCTION
    arr_t<arr_t<unsigned int, Trigon<double, P2, Dim>::count>,
               Trigon<double, P1, Dim>::count>
    calculate_f()
    {
        //scoped_simple_timer("trigon_cal_f");
        arr_t<arr_t<unsigned int, Trigon<double, P2, Dim>::count>,
                   Trigon<double, P1, Dim>::count>
            retval;

        for (unsigned int i = 0; i < Trigon<double, P1, Dim>::count; ++i) {
            for (unsigned int j = 0; j < Trigon<double, P2, Dim>::count; ++j) {

#if 0
                auto indices1(indices<P1, Dim>());
                Index_t<P1> index1(indices1[i]);
                auto indices2(indices<P2, Dim>());
                Index_t<P2> index2(indices2[j]);
#endif
#if 0
                Index_t<P1> index1 = canonical_to_index<P1, Dim>()[i];
                Index_t<P2> index2 = canonical_to_index<P2, Dim>()[j];
#endif
#if 1
                Index_t<P1> index1 = indices<P1, Dim>()[i];
                Index_t<P2> index2 = indices<P2, Dim>()[j];
#endif

                Index_t<P1 + P2> index3;

                size_t m = 0;
                for (; m < P1; ++m) {
                    index3[m] = index1[m];
                }
                for (size_t n = 0; n < P2; ++m, ++n) {
                    index3[m] = index2[n];
                }
                std::sort(index3.begin(), index3.end());
                //trigon_impl::sort(index3);
                retval.at(i).at(j) = index_to_canonical<P1 + P2, Dim>().at(index3);
            }
        }
        return retval;
    }

#if 0
    template <unsigned int P1, unsigned int P2>
    KOKKOS_INLINE_FUNCTION
    unsigned int f(unsigned int i, unsigned j)
    {
#if 1
        static arr_t<
            arr_t<unsigned int, Trigon<double, P2, Dim>::count>,
            Trigon<double, P1, Dim>::count>
            mapping = calculate_f<P1, P2>();
        return mapping[i][j];
#endif

#if 0
        auto mapping = calculate_f<P1, P2>();
        return mapping[i][j];
#endif
    }
#endif

    template <unsigned int P2>
    KOKKOS_INLINE_FUNCTION
    unsigned int f(unsigned int i, unsigned j)
    {
#ifdef __CUDA_ARCH__
        return 0;
#else
        static arr_t<
            arr_t<unsigned int, Trigon<double, P2, Dim>::count>,
            Trigon<double, Power, Dim>::count>
            mapping = calculate_f<Power, P2>();
        return mapping[i][j];
#endif
    }


    template <unsigned int New_power, typename Mult_trigon_t, typename Array_t>
    KOKKOS_INLINE_FUNCTION
    void collect_products(Mult_trigon_t const& t, Array_t& new_terms)
    {
        //simple_timer_start("trigon_collect_products");
        // this x right = new
        constexpr unsigned int right_power = New_power - Power;
        if (right_power <= t.power()) {
            auto& right_terms(t.template get_subpower<right_power>().terms);
            for (size_t i = 0; i < terms.size(); ++i) {
                for (size_t j = 0; j < right_terms.size(); ++j) {
                    //size_t k = f<Power, right_power>(i, j);
                    size_t k = f<right_power>(i, j);
                    new_terms[k] += terms[i] * right_terms[j];
                }
            }
        }
        //simple_timer_stop("trigon_collect_products");
        lower.template collect_products<New_power>(t, new_terms);
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator*=(Trigon<T, Power, Dim> const& t)
    {
        if (Power > 1) {
            //simple_timer_start("trigon_*=(T)");
            Terms_t new_terms;
            new_terms.fill(0);
            collect_products<Power>(t, new_terms);
            //simple_timer_stop("trigon_*=(T)");
            lower *= t.lower;
            terms = new_terms;
        } else {
            //scoped_simple_timer("trigon_*=(T)");
            const T this_value = value();
            const T right_value = t.value();
            for (size_t i = 0; i < Dim; ++i) {
                terms[i] *= right_value;
                terms[i] += t.terms[i] * this_value;
            }
            value() *= right_value;
        }

        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator*=(T val)
    {
        //simple_timer_start("trigon_*=(d)");
        for (auto&& c : terms) {
            c *= val;
        }
        //simple_timer_stop("trigon_*=(d)");
        lower *= val;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator*(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval *= t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator*(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval *= val;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator/=(Trigon<T, Power, Dim> const& t)
    {
        // this / t = new
        if (Power > 1) {
            lower /= t.lower;
            Terms_t new_terms;
            new_terms.fill(0);
            lower.template collect_products<Power>(t, new_terms);
            T t0 = t.value();
            for (size_t i = 0; i < new_terms.size(); ++i) {
                terms[i] = (terms[i] - new_terms[i]) / t0;
            }
        } else {
            const T this_value = value();
            const T right_value = t.value();
            const T inv_right_value2 = 1.0 / (right_value * right_value);
            for (size_t i = 0; i < Dim; ++i) {
                terms[i] *= right_value;
                terms[i] -= t.terms[i] * this_value;
                terms[i] *= inv_right_value2;
            }
            value() /= right_value;
        }
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator/=(T val)
    {
        for (auto&& c : terms) {
            c /= val;
        }
        lower /= val;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator/(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval /= t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, Power, Dim> operator/(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval /= val;
        return retval;
    }

    // partial derivative
    // [0, 0] => dTrigon/(dx dx)
    // [0, 1] => dTrigon/(dx dy)
    // [2, 2, 2] => dTrigon/(dz dz dz)
    template<size_t DP>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if_t<(Power>DP && DP>1), Trigon<T, Power-DP, Dim>>
    derivative(arr_t<unsigned int, DP> const& idx)
    {
        arr_t<unsigned int, DP-1> i2;
        for(int i=0; i<DP-1; ++i) i2[i] = idx[i];
        return partial_deriv(*this, idx[DP-1]).derivative(i2);
    }

    template<size_t DP>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if_t<(Power>DP && DP==1), Trigon<T, Power-DP, Dim>>
    derivative(arr_t<unsigned int, DP> const& idx)
    {
        return partial_deriv(*this, idx[0]);
    }


    // evaluation
    KOKKOS_INLINE_FUNCTION
    T operator()(arr_t<T, Dim> const& x) const
    {
        T val{};

        // current power
        auto inds = indices<Power, Dim>();
        for(int i=0; i<terms.size(); ++i)
        {
            if (abs(terms[i]))
            {
                T t = terms[i];
                for(auto idx : inds[i]) t *= x[idx];
                val += t;
            }
        }

        // lower power
        val += lower(x);

        return val;
    }

    // composition
    template<unsigned int P>
    KOKKOS_INLINE_FUNCTION
    Trigon<T, P, Dim> compose(TMapping<Trigon<T, P, Dim>> const& x) const
    {
        Trigon<T, P, Dim> val;

        // current power
        auto inds = indices<Power, Dim>();
        for(int i=0; i<terms.size(); ++i)
        {
            if (abs(terms[i]))
            {
                Trigon<T, P, Dim> t(terms[i]);
                for(auto idx : inds[i]) t *= x[idx];
                val += t;
            }
        }

        // lower power
        val += lower.compose(x);

        return val;
    };

#ifdef __CUDA_ARCH__

    syn::dummy_json to_json() const
    { return syn::dummy_json{}; }

#else

    syn::json to_json() const
    { 
        syn::json val = {
            { "dim", Dim },
            { "power", Power },
            { "terms", syn::json::array() }
        }; 

        syn::json terms{};
        to_json_impl(terms); 

        val["terms"] = std::move(terms);

        return val; 
    }

    void to_json_impl(syn::json& val) const
    {
        // lower power
        lower.to_json_impl(val);

        // current power
        syn::json v{};
        v["power"] = Power;
        v["terms"] = syn::json::array();

        auto inds = indices<Power, Dim>();
        for(int i=0; i<terms.size(); ++i)
        {
            arr_t<int, Dim> exp;
            for(auto idx : inds[i]) ++exp[idx];

            syn::json term = {
                {"exp", syn::json::array()},
                {"term", term_to_json_val(terms[i])}
            };

            for(int i=0; i<exp.size(); ++i) 
                term["exp"][i] = exp[i];

            v["terms"].emplace_back(std::move(term));
        }

        val.emplace_back(std::move(v));
    }

#endif

    friend std::ostream& operator<<(std::ostream& os, 
            Trigon<T, Power, Dim> const& t)
    {
        os << t.lower << "P(" << Power << "): (";
        //for(int i=0; i<t.count-1; ++i) os << std::setprecision(16) << t.terms[i] << ", ";
        for(int i=0; i<t.count-1; ++i) os << t.terms[i] << ", ";
        os << t.terms[t.count-1] << ")\n";
        return os;
    }

    template<class AR>
    void serialize(AR& ar)
    {
        ar(lower);
        ar(terms);
    }
};

template<typename T, unsigned int P, unsigned int D>
struct is_trigon<Trigon<T, P, D>> : std::true_type { };

// power 0
template <typename T, unsigned int Dim>
class Trigon<T, 0, Dim>
{
public:
    static constexpr unsigned int count = 1;
    typedef arr_t<T, count> Terms_t;
    Terms_t terms;

    // implicit conversion
    template<class U>
    KOKKOS_INLINE_FUNCTION
    explicit Trigon(Trigon<U, 0, Dim> const& o)
    { terms.from(o.terms); }

    KOKKOS_INLINE_FUNCTION
    Trigon() { terms[0] = 0; }

    KOKKOS_INLINE_FUNCTION
    Trigon(T val) { terms[0] = val; }

    KOKKOS_INLINE_FUNCTION
    static constexpr unsigned int power() { return 0; }

    KOKKOS_INLINE_FUNCTION
    const T value() const { return terms[0]; }

    KOKKOS_INLINE_FUNCTION
    void set(T val) { terms[0] = val; }

    template <unsigned int Subpower>
    KOKKOS_INLINE_FUNCTION
    Trigon<T, Subpower, Dim>& get_subpower()
    { return *this; }

    template <unsigned int Subpower>
    KOKKOS_INLINE_FUNCTION
    const Trigon<T, Subpower, Dim>& get_subpower() const
    { return *this; }

    KOKKOS_INLINE_FUNCTION
    bool operator== (T rhs) const
    { 
        double const eps = 1e5 * std::numeric_limits<double>::epsilon();
        return !(abs(terms[0] - rhs) > eps);
    }

    KOKKOS_INLINE_FUNCTION
    bool operator!= (T rhs) const
    { return !((*this) == rhs); }

    template <typename F>
    KOKKOS_INLINE_FUNCTION
    void each_term(F f)
    { f(0, arr_t<size_t, 0>(), terms[0]); }

    KOKKOS_INLINE_FUNCTION
    void set_term(size_t idx, T const& val)
    { terms[0] = val; }

    KOKKOS_INLINE_FUNCTION
    T get_term(size_t idx)
    { return terms[0]; }

    KOKKOS_INLINE_FUNCTION
    void set_term(unsigned int power, size_t idx, T const& val)
    { set_term(idx, val); }

    KOKKOS_INLINE_FUNCTION
    T get_term(unsigned int power, size_t idx)
    { return get_term(idx); }

    // keep terms with power in [lower, upper]
    KOKKOS_INLINE_FUNCTION
    void filter(unsigned int p_lower, unsigned int p_upper)
    {
        if (p_lower==0) terms[0] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim>& operator+=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] += t.terms[0];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim>& operator+=(T t)
    {
        terms[0] += t;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator+(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval += t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator+(T t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval += t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator-() const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval.terms[0] = -terms[0];
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim>& operator-=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] -= t.terms[0];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim>& operator-=(T t)
    {
        terms[0] -= t;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator-(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval -= t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator-(T t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval -= t;
        return retval;
    }

    template <unsigned int New_power, typename Mult_trigon_t, typename Array_t>
    KOKKOS_INLINE_FUNCTION
    void collect_products(Mult_trigon_t const& t, Array_t& new_terms)
    {
        // this x right = new
        constexpr unsigned int right_power = New_power;
        if (right_power <= t.power()) {
            auto& right_terms(t.template get_subpower<right_power>().terms);
            for (size_t j = 0; j < right_terms.size(); ++j) {
                new_terms[j] += terms[0] * right_terms[j];
            }
        }
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator*=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] *= t.terms[0];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator*=(T val)
    {
        terms[0] *= val;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator*(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval *= t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator*(T val) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval *= val;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator/=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] /= t.terms[0];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator/=(T val)
    {
        terms[0] /= val;
        return *this;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator/(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval /= t;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    Trigon<T, 0, Dim> operator/(T val) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval /= val;
        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    T operator()(arr_t<T, Dim> const& x) const
    { return terms[0]; }

    template<unsigned int P>
    KOKKOS_INLINE_FUNCTION
    Trigon<T, P, Dim> compose(TMapping<Trigon<T, P, Dim>> const& x) const
    { return Trigon<T, P, Dim>(terms[0]); };

    friend std::ostream& 
    operator<<(std::ostream& os, Trigon<T, 0, Dim> const& t)
    { os << "\nP(0): (" << t.terms[0] << ")\n"; return os; }

#ifdef __CUDA_ARCH__

    syn::dummy_json to_json() const
    { return syn::dummy_json{}; }

#else

    syn::json to_json() const
    { 
        syn::json val = {
            { "dim", Dim },
            { "power", 0},
            { "terms", syn::json::array() }
        }; 

        syn::json terms{};
        to_json_impl(terms); 

        val["terms"] = std::move(terms);

        return val; 
    }

    void to_json_impl(syn::json& val) const
    {
        // current power
        syn::json v{};
        v["power"] = 0;
        v["terms"][0] = {
            {"exp", syn::json::array()},
            {"term", term_to_json_val(terms[0])}
        };

        for(int i=0; i<Dim; ++i)
            v["terms"][0]["exp"][i] = 0;

        val.emplace_back(std::move(v));
    }

#endif

    template<class AR>
    void serialize(AR& ar)
    {
        ar(terms);
    }
};


template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
operator+(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval(t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, Power, Dim>
operator+(Trigon<T, Power, Dim> const& t1,
          Trigon<std::complex<T>, Power, Dim> const& t2)
{
    Trigon<std::complex<T>, Power, Dim> retval;
    for (size_t i = 0; i < retval.terms.size(); ++i) {
        retval.terms[i] = t1.terms[i] + t2.terms[i];
    }
    retval.lower = t1.lower + t2.lower;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, Power, Dim>
operator+(Trigon<std::complex<T>, Power, Dim> const& t1,
          Trigon<T, Power, Dim> const& t2)
{
    return t2 + t1;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, 0, Dim>
operator+(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval(t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, 0, Dim> operator+(
    Trigon<T, 0, Dim> const& t1, Trigon<std::complex<T>, 0, Dim> const& t2)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = t1.terms[0] + t2.terms[0];
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
operator-(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval(-t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, 0, Dim>
operator-(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval(-t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, Power, Dim>
operator-(Trigon<T, Power, Dim> const& t1,
          Trigon<std::complex<T>, Power, Dim> const& t2)
{
    Trigon<std::complex<T>, Power, Dim> retval;
    for (size_t i = 0; i < retval.terms.size(); ++i) {
        retval.terms[i] = t1.terms[i] - t2.terms[i];
    }
    retval.lower = t1.lower - t2.lower;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, Power, Dim>
operator-(Trigon<std::complex<T>, Power, Dim> const& t1,
          Trigon<T, Power, Dim> const& t2)
{
    Trigon<std::complex<T>, Power, Dim> retval;
    for (size_t i = 0; i < retval.terms.size(); ++i) {
        retval.terms[i] = t1.terms[i] - t2.terms[i];
    }
    retval.lower = t1.lower - t2.lower;
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, 0, Dim> operator-(
    Trigon<T, 0, Dim> const& t1, Trigon<std::complex<T>, 0, Dim> const& t2)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = t1.terms[0] - t2.terms[0];
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, 0, Dim> operator-(
    Trigon<std::complex<T>, 0, Dim> const& t1, Trigon<T, 0, Dim> const& t2)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = t1.terms[0] - t2.terms[0];
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim> operator*(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval(t);
    retval *= val;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, Power, Dim> 
operator*(std::complex<T> val, Trigon<T, Power, Dim> const& t)
{
    Trigon<std::complex<T>, Power, Dim> retval;
    for (size_t i = 0; i < retval.terms.size(); ++i) {
        retval.terms[i] = val * t.terms[i];
    }
    retval.lower = val * t.lower;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, Power, Dim> 
operator*( Trigon<std::complex<T>, Power, Dim> const& t1,
    Trigon<T, Power, Dim> const& t2)
{
    return t1 * (std::complex<T>(1.0, 0.0) * t2);
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, Power, Dim> 
operator*( Trigon<T, Power, Dim> const& t1,
    Trigon<std::complex<T>, Power, Dim> const& t2)
{
    return (std::complex<T>(1.0, 0.0) * t1) * t2;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, 0, Dim> operator*(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval(t);
    retval *= val;
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<std::complex<T>, 0, Dim> operator*(std::complex<T> val,
                                          Trigon<T, 0, Dim> const& t)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = val * t.terms[0];
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
operator/(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval;
    retval.template get_subpower<0>().terms[0] = val;
    retval /= t;
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, 0, Dim>
operator/(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval;
    retval.terms[0] = val / t.terms[0];
    return retval;
}

#if 1
template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
typename std::enable_if<(Power == 1), void>::type
calculate_partial(Trigon<T, Power, Dim> const& source, size_t index,
                  Trigon<T, Power - 1, Dim>& dest)
{
    dest.terms[0] = source.terms[index];
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
typename std::enable_if<(Power > 1), void>::type
calculate_partial(Trigon<T, Power, Dim> const& source, size_t index,
                  Trigon<T, Power - 1, Dim>& dest)
{
    auto source_indices = indices<Power, Dim>();

    for (size_t i = 0; i < source_indices.size(); ++i) 
    {
        auto& source_index = source_indices[i];
        Index_t<Power - 1> dest_index;
        size_t exponent = 0;
        size_t k = 0;

        for (size_t j = 0; j < source_index.size(); ++j) {
            if (source_index[j] == index) {
                exponent += 1;
            } else {
                dest_index[k] = source_index[j];
                ++k;
            }
        }

        if (exponent > 0) {
            for (size_t j = 0; j < (exponent - 1); ++j, ++k) {
                dest_index[k] = index;
            }
            std::sort(dest_index.begin(), dest_index.end());
            //trigon_impl::sort(dest_index);
            dest.terms[index_to_canonical<Power - 1, Dim>().at(dest_index)] =
                T(exponent) * source.terms[i];
        }
    }

    calculate_partial<T, Power - 1, Dim>(source.lower, index, dest.lower);
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power - 1, Dim>
partial_deriv(Trigon<T, Power, Dim> const& t, size_t index)
{
    Trigon<T, Power - 1, Dim> retval;
    calculate_partial<T, Power, Dim>(t, index, retval);
    return retval;
}
#endif

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
real(Trigon<std::complex<T>, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval;
    for (size_t i = 0; i < retval.terms.size(); ++i) {
        retval.terms[i] = t.terms[i].real();
    }
    retval.lower = real(t.lower);
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, 0, Dim> 
real(Trigon<std::complex<T>, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval;
    retval.terms[0] = t.terms[0].real();
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
imag(Trigon<std::complex<T>, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval;
    for (size_t i = 0; i < retval.terms.size(); ++i) {
        retval.terms[i] = t.terms[i].imag();
    }
    retval.lower = imag(t.lower);
    return retval;
}

template <typename T, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, 0, Dim> 
imag(Trigon<std::complex<T>, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval;
    retval.terms[0] = t.terms[0].imag();
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
exp(Trigon<T, Power, Dim> const& t)
{
    T val = t.value();
    Trigon<T, Power, Dim> x(t - val);
    Trigon<T, Power, Dim> xn(1);
    Trigon<T, Power, Dim> retval(1.0);
    double fact_n = 1;
    for (size_t n = 1; n <= Power; ++n) {
        xn *= x;
        fact_n *= n;
        retval += xn / fact_n;
    }
    retval *= std::exp(val);
    return retval;
}

constexpr int max_power = 7;

template <typename T>
using Derivatives_t = arr_t<T, max_power + 1>;

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
generic_transcendental(Trigon<T, Power, Dim> const& t,
                       Derivatives_t<T> const& derivatives)
{
    T val = t.value();
    Trigon<T, Power, Dim> x(t - val);
    Trigon<T, Power, Dim> xn(1);
    Trigon<T, Power, Dim> retval(1.0);
    int fact_n = 1;
    for (size_t n = 1; n <= Power; ++n) {
        xn *= x;
        fact_n *= n;
        retval += xn * derivatives.at(n) / fact_n;
    }
    retval.value() = derivatives[0];
    return retval;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
Derivatives_t<T>
sin_derivatives(T x, unsigned int power)
{
    Derivatives_t<T> retval;
    retval[0] = std::sin(x);
    if (power > 0) {
        retval[1] = std::cos(x);
    }
    if (power > 1) {
        retval[2] = -retval[0];
    }
    if (power > 2) {
        retval[3] = -retval[1];
    }
    if (power > 3) {
        retval[4] = retval[0];
    }
    if (power > 4) {
        retval[5] = retval[1];
    }
    if (power > 5) {
        retval[6] = retval[2];
    }
    if (power > 6) {
        retval[7] = retval[3];
    }
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
sin(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, sin_derivatives(t.value(), Power));
}

template <typename T>
KOKKOS_INLINE_FUNCTION
Derivatives_t<T>
cos_derivatives(T x, unsigned int power)
{
    Derivatives_t<T> retval;
    retval[0] = std::cos(x);
    if (power > 0) {
        retval[1] = -std::sin(x);
    }
    if (power > 1) {
        retval[2] = -retval[0];
    }
    if (power > 2) {
        retval[3] = -retval[1];
    }
    if (power > 3) {
        retval[4] = retval[0];
    }
    if (power > 4) {
        retval[5] = retval[1];
    }
    if (power > 5) {
        retval[6] = retval[2];
    }
    if (power > 6) {
        retval[7] = retval[3];
    }
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
cos(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, cos_derivatives(t.value(), Power));
}

template <typename T>
KOKKOS_INLINE_FUNCTION
Derivatives_t<T>
tan_derivatives(T x, unsigned int power)
{
    // for qpow
    using namespace kt;

    Derivatives_t<T> retval;
    T tanx(std::tan(x));
    T secx;
    retval[0] = tanx;
    if (power > 0) {
        secx = 1 / std::cos(x);
        retval[1] = qpow(secx, 2);
    }
    if (power > 1) {
        retval[2] = 2 * tanx * qpow(secx, 2);
    }
    if (power > 2) {
        retval[3] = 4 * qpow(secx, 2) * qpow(tanx, 2) + 2 * qpow(secx, 4);
    }
    if (power > 3) {
        retval[4] =
            8 * qpow(secx, 2) * qpow(tanx, 3) + 16 * tanx * qpow(secx, 4);
    }
    if (power > 4) {
        retval[5] = 16 * qpow(secx, 2) * qpow(tanx, 4) +
                    88 * qpow(secx, 4) * qpow(tanx, 2) + 16 * qpow(secx, 6);
    }
    if (power > 5) {
        retval[6] = 32 * qpow(secx, 2) * qpow(tanx, 5) +
                    416 * qpow(secx, 4) * qpow(tanx, 3) +
                    272 * tanx * qpow(secx, 6);
    }
    if (power > 6) {
        retval[7] = 64 * qpow(secx, 2) * qpow(tanx, 6) +
                    1824 * qpow(secx, 4) * qpow(tanx, 4) +
                    2880 * qpow(secx, 6) * qpow(tanx, 2) + 272 * qpow(secx, 8);
    }
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<T, Power, Dim>
tan(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, tan_derivatives(t.value(), Power));
}

template <typename T>
KOKKOS_INLINE_FUNCTION
Derivatives_t<T>
sqrt_derivatives(T x, unsigned int power)
{
    Derivatives_t<T> retval;
    T invx;
    retval[0] = std::sqrt(x);
    if (power > 0) {
        retval[1] = 0.5 / retval[0];
    }
    if (power > 1) {
        invx = 1 / x;
        retval[2] = -0.5 * invx * retval[1];
    }
    if (power > 2) {
        retval[3] = -1.5 * invx * retval[2];
    }
    if (power > 3) {
        retval[4] = -2.5 * invx * retval[3];
    }
    if (power > 4) {
        retval[5] = -3.5 * invx * retval[4];
    }
    if (power > 5) {
        retval[6] = -4.5 * invx * retval[5];
    }
    if (power > 6) {
        retval[7] = -5.5 * invx * retval[6];
    }
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
sqrt(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, sqrt_derivatives(t.value(), Power));
}

template <typename T>
KOKKOS_INLINE_FUNCTION
Derivatives_t<T>
asin_derivatives(T x, unsigned int power)
{
    // qpow
    using namespace kt;

    Derivatives_t<T> retval;
    retval[0] = std::asin(x);
    T x2(x * x);
    T invsqrt1mx2(1.0 / std::sqrt(1 - x2));
    if (power > 0) {
        retval[1] = invsqrt1mx2;
    }
    if (power > 1) {
        retval[2] = qpow(invsqrt1mx2, 3) * x;
    }
    if (power > 2) {
        retval[3] = 2 * qpow(invsqrt1mx2, 5) * x2 + qpow(invsqrt1mx2, 5);
    }
    if (power > 3) {
        retval[4] =
            6 * qpow(invsqrt1mx2, 7) * x * x2 + 9 * qpow(invsqrt1mx2, 7) * x;
    }
    if (power > 4) {
        retval[5] = 24 * qpow(invsqrt1mx2, 9) * qpow(x2, 2) +
                    72 * qpow(invsqrt1mx2, 9) * x2 + 9 * qpow(invsqrt1mx2, 9);
    }
    if (power > 5) {
        retval[6] = 120 * qpow(invsqrt1mx2, 11) * x * qpow(x2, 2) +
                    600 * qpow(invsqrt1mx2, 11) * x * x2 +
                    225 * qpow(invsqrt1mx2, 11) * x;
    }
    if (power > 6) {
        retval[7] = 720 * qpow(invsqrt1mx2, 13) * qpow(x2, 3) +
                    5400 * qpow(invsqrt1mx2, 13) * qpow(x2, 2) +
                    4050 * qpow(invsqrt1mx2, 13) * x2 +
                    225 * qpow(invsqrt1mx2, 13);
    }
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
asin(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, asin_derivatives(t.value(), Power));
}

template <typename T>
KOKKOS_INLINE_FUNCTION
Derivatives_t<T>
atan_derivatives(T x, unsigned int power)
{
    //throw std::runtime_error("atan(trigon t) yet to be implemented");
    return Derivatives_t<T>();

#if 0
    Derivatives_t<T> retval;
    retval[0] = std::atan(x);
    T x2(x * x);
    T invsqrt1mx2(1.0 / std::sqrt(1 - x2));
    if (power > 0) {
        retval[1] = invsqrt1mx2;
    }
    if (power > 1) {
        retval[2] = qpow(invsqrt1mx2, 3) * x;
    }
    if (power > 2) {
        retval[3] = 2 * qpow(invsqrt1mx2, 5) * x2 + qpow(invsqrt1mx2, 5);
    }
    if (power > 3) {
        retval[4] =
            6 * qpow(invsqrt1mx2, 7) * x * x2 + 9 * qpow(invsqrt1mx2, 7) * x;
    }
    if (power > 4) {
        retval[5] = 24 * qpow(invsqrt1mx2, 9) * qpow(x2, 2) +
                    72 * qpow(invsqrt1mx2, 9) * x2 + 9 * qpow(invsqrt1mx2, 9);
    }
    if (power > 5) {
        retval[6] = 120 * qpow(invsqrt1mx2, 11) * x * qpow(x2, 2) +
                    600 * qpow(invsqrt1mx2, 11) * x * x2 +
                    225 * qpow(invsqrt1mx2, 11) * x;
    }
    if (power > 6) {
        retval[7] = 720 * qpow(invsqrt1mx2, 13) * qpow(x2, 3) +
                    5400 * qpow(invsqrt1mx2, 13) * qpow(x2, 2) +
                    4050 * qpow(invsqrt1mx2, 13) * x2 +
                    225 * qpow(invsqrt1mx2, 13);
    }
    return retval;
#endif
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
atan(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, atan_derivatives(t.value(), Power));
}


template <typename T>
KOKKOS_INLINE_FUNCTION
Derivatives_t<T>
log_derivatives(T x, unsigned int power)
{
    Derivatives_t<T> retval;
    retval[0] = std::log(x);
    T invx(1.0 / x);
    T powinvx = 1.0;
    T fact = 1.0;
    for (size_t i = 1; i < power; ++i) {
        powinvx *= invx;
        if (i > 1) {
            fact *= -(i - 1);
        }
        retval[i] = fact * powinvx;
    }
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
KOKKOS_INLINE_FUNCTION
Trigon<T, Power, Dim>
log(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, log_derivatives(t.value(), Power));
}



template<typename TRIGON>
struct TMapping
{
    using trigon_t = TRIGON;
    constexpr static unsigned int dim = TRIGON::dim;
    constexpr static unsigned int power = TRIGON::power();

    //template<size_t ROW, size_t OPT>
    using matrix_t = Eigen::Matrix<typename trigon_t::data_type, Eigen::Dynamic, dim>;

    arr_t<TRIGON, TRIGON::dim> comp;

    template<class U>
    KOKKOS_INLINE_FUNCTION
    explicit TMapping(TMapping<Trigon<U, power, dim>> const& o)
    { comp.from(o.comp); }

    KOKKOS_INLINE_FUNCTION
    TMapping() : comp()
    { }

    KOKKOS_INLINE_FUNCTION
    TRIGON& operator[](size_t idx) { return comp[idx]; }

    KOKKOS_INLINE_FUNCTION
    TRIGON const& operator[](size_t idx) const { return comp[idx]; }

    // evaluation
    KOKKOS_INLINE_FUNCTION
    arr_t<typename TRIGON::data_type, dim>
    operator()(arr_t<typename TRIGON::data_type, dim> const& x) const
    {
        arr_t<typename TRIGON::data_type, dim> val;
        for(int i=0; i<dim; ++i) val[i] = comp[i](x);
        return val;
    }

    // compose
    KOKKOS_INLINE_FUNCTION
    TMapping<TRIGON>
    operator()(TMapping<TRIGON> const& x, 
            arr_t<typename TRIGON::data_type, dim> const& ref = {0}) const
    {
        TMapping<TRIGON> ret;
        TMapping<TRIGON> u = x;

        for(int i=0; i<comp.size(); ++i) 
            u[i] = x[i] - TRIGON(ref[i]);

        for(int i=0; i<comp.size(); ++i) 
            ret[i] = comp[i].compose(u);

        return ret;
    }

    // subtract
    KOKKOS_INLINE_FUNCTION
    TMapping<TRIGON>
    operator-(TMapping<TRIGON> const& x)
    {
        TMapping<TRIGON> ret;
        for(int i=0; i<dim; ++i) ret[i] = comp[i] - x.comp[i];
        return ret;
    }

    // exp
    KOKKOS_INLINE_FUNCTION
    TRIGON operator^(TRIGON const& x)
    {
        TRIGON answer;
        const int s = x.dim; // space dim

        for(int i=0; i<s; ++i)
        {
            //answer += comp[i] * partial_deriv(x, i);

            auto xd = partial_deriv(x, i);
            TRIGON xdp; xdp.lower = xd;

            answer += comp[i] * xdp;
        }

        return answer;
    }

    KOKKOS_INLINE_FUNCTION
    TRIGON
    exp_map(typename TRIGON::data_type const& t, TRIGON const& x)
    {
        constexpr int MX_MAXITER = 100;
        typename TRIGON::data_type zero{};

        double f = 1.0;
        int count = 0;

        TRIGON answer = x;
        TRIGON u = (t/f) * ( (*this)^x );
        ++f;

        while( ++count<MX_MAXITER && u!=zero )
        {
            answer += u;
            u = (t/f) * ( (*this)^u );
            ++f;
        }

        if (count >= MX_MAXITER)
        {
            std::cout << "TMapping::exp_map() number of iteration has "
                "exceeded " << MX_MAXITER << " without achieving "
                "convergence. Results maybe incorrect.\n";
        }

        return answer;
    }


    KOKKOS_INLINE_FUNCTION
    TMapping<TRIGON>
    exp_map(typename TRIGON::data_type const& t, TMapping<TRIGON> const& m)
    {
        TMapping<TRIGON> z;
        for(int i=0; i<dim; ++i) z[i] = exp_map(t, m[i]);
        return z;
    }


    KOKKOS_INLINE_FUNCTION
    void filter(unsigned int lower, unsigned int upper)
    {
        for(auto & t : comp) t.filter(lower, upper);
    }

    KOKKOS_INLINE_FUNCTION
    karray2d_row jacobian() const
    {
        karray2d_row ret("jacobian", dim, dim);
        for(int i=0; i<dim; ++i)
            for(int j=0; j<dim; ++j)
                ret(i,j) = comp[i].template get_subpower<1>().terms[j];
        return ret;
    }

#ifdef __CUDA_ARCH__
    syn::dummy_json to_json() const
    { return syn::dummy_json{}; }
#else
    syn::json to_json() const
    { 
        syn::json val = syn::json::array();
        for(auto const& t : comp) 
            val.emplace_back(t.to_json());
        return val; 
    }
#endif

    template<class AR>
    void serialize(AR& ar)
    {
        ar(comp);
    }

};

#if 1
template<class TRIGON>
std::ostream& operator<<(
        std::ostream& os, TMapping<TRIGON> const& m)
{
    for(auto const& t : m.comp) os << t;
    return os;
}


template<typename TRIGON>
TMapping<TRIGON> operator*(
        typename TMapping<TRIGON>::matrix_t const& m,
        TMapping<TRIGON> const& x)
{
    TMapping<TRIGON> z;

    for(int i=0; i< m.rows(); ++i)
    {
        z[i] = m(i,0) * x[0];

        int j = 1;
        while(j < m.cols())
        {
            z[i] += m(i,j) * x[j];
            ++j;
        }
    }

    return z;
}
#endif



#endif // TRIGON_H
