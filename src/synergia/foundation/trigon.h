#ifndef TRIGON_H
#define TRIGON_H

#include <array>
#include <complex>
#include <iostream> // jfa remove me!!!!
#include <type_traits>
#include <unordered_map>

constexpr int
factorial(int n)
{
    return n <= 1 ? 1 : (n * factorial(n - 1));
}

template <typename T, unsigned int Power, unsigned int Dim>
class Trigon;

template <unsigned int Power>
using Index_t = std::array<size_t, Power>;

template <unsigned int Power, unsigned int Dim>
using Indices_t = std::array<Index_t<Power>, Trigon<double, Power, Dim>::count>;

template <unsigned int Power, unsigned int Dim>
typename std::enable_if<((Power == 1) || (Power == 0)),
                        Indices_t<Power, Dim>>::type
indices()
{
    Indices_t<Power, Dim> retval;
    for (size_t i = 0; i < Dim; ++i) {
        retval[0][i] = i;
    }

    return retval;
}

template <unsigned int Power, unsigned int Dim>
typename std::enable_if<(Power > 1), Indices_t<Power, Dim>>::type
indices()
{
    Indices_t<Power - 1, Dim> subindices = indices<Power - 1, Dim>();
    Indices_t<Power, Dim> retval;
    size_t count = 0;
    for (size_t i = 0; i < Dim; ++i) {
        for (size_t j = 0; j < subindices.size(); ++j) {
            if (subindices[j][0] >= i) {
                std::array<size_t, Power> entry;
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
    size_t operator()(std::array<size_t, Length> const& arr) const
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
using Map_t =
    std::unordered_map<std::array<size_t, Power>, size_t, Array_hash<Power>>;

template <unsigned int Power, unsigned int Dim>
const Map_t<Power, Dim>
fill_index_to_canonical()
{
    Map_t<Power, Dim> map;
    Indices_t<Power, Dim> the_indices = indices<Power, Dim>();
    for (size_t i = 0; i < the_indices.size(); ++i) {
        map[the_indices[i]] = i;
    }

    return map;
}

template <unsigned int Power, unsigned int Dim>
const Map_t<Power, Dim>&
index_to_canonical()
{
    static const Map_t<Power, Dim> map = fill_index_to_canonical<Power, Dim>();
    return map;
}

template <typename T, unsigned int Power, unsigned int Dim>
class Trigon
{
public:
    Trigon<T, Power - 1, Dim> lower;
    static constexpr unsigned int count =
        factorial(Dim + Power - 1) / (factorial(Dim - 1) * factorial(Power));
    typedef std::array<T, count> Terms_t;
    Terms_t terms;

    Trigon() { terms.fill(0); }

    Trigon(T val)
        : lower(val)
    {
        terms.fill(0);
    }

    Trigon(T val, size_t index)
        : lower(val)
    {
        terms.fill(0);
        get_subpower<1>().terms[index] = 1; // jfa fixme
    }

    static constexpr unsigned int power() { return Power; }

    const T& value() const { return get_subpower<0>().terms[0]; }

    T& value() { return get_subpower<0>().terms[0]; }

    template <unsigned int Subpower>
    typename std::enable_if<(Subpower < Power), Trigon<T, Subpower, Dim>&>::type
    get_subpower()
    {
        return lower.template get_subpower<Subpower>();
    }

    template <unsigned int Subpower>
    const typename std::enable_if<(Subpower < Power),
                                  Trigon<T, Subpower, Dim> const&>::type
    get_subpower() const
    {
        return lower.template get_subpower<Subpower>();
    }

    template <unsigned int Subpower>
    typename std::enable_if<(Subpower == Power),
                            Trigon<T, Subpower, Dim>&>::type
    get_subpower()
    {
        return *this;
    }

    template <unsigned int Subpower>
    const typename std::enable_if<(Subpower == Power),
                                  Trigon<T, Subpower, Dim> const&>::type
    get_subpower() const
    {
        return *this;
    }

    Trigon<T, Power, Dim>& operator+=(Trigon<T, Power, Dim> const& t)
    {
        lower += t.lower;
        for (size_t i = 0; i < terms.size(); ++i) {
            terms[i] += t.terms[i];
        }
        return *this;
    }

    Trigon<T, Power, Dim>& operator+=(T val)
    {
        lower += val;
        return *this;
    }

    Trigon<T, Power, Dim> operator+(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval += t;
        return retval;
    }

    Trigon<T, Power, Dim> operator+(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval += val;
        return retval;
    }

    Trigon<T, Power, Dim> operator-() const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval.lower = -lower;
        for (size_t i = 0; i < terms.size(); ++i) {
            retval.terms[i] = -terms[i];
        }
        return retval;
    }

    Trigon<T, Power, Dim>& operator-=(Trigon<T, Power, Dim> const& t)
    {
        lower -= t.lower;
        for (size_t i = 0; i < terms.size(); ++i) {
            terms[i] -= t.terms[i];
        }
        return *this;
    }

    Trigon<T, Power, Dim>& operator-=(T val)
    {
        lower -= val;
        return *this;
    }

    Trigon<T, Power, Dim> operator-(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval -= t;
        return retval;
    }

    Trigon<T, Power, Dim> operator-(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval -= val;
        return retval;
    }

    template <unsigned int P1, unsigned int P2>
    std::array<std::array<unsigned int, Trigon<double, P2, Dim>::count>,
               Trigon<double, P1, Dim>::count>
    calculate_f()
    {
        std::array<std::array<unsigned int, Trigon<double, P2, Dim>::count>,
                   Trigon<double, P1, Dim>::count>
            retval;

        for (unsigned int i = 0; i < Trigon<double, P1, Dim>::count; ++i) {
            for (unsigned int j = 0; j < Trigon<double, P2, Dim>::count; ++j) {

                auto indices1(indices<P1, Dim>());
                Index_t<P1> index1(indices1[i]);
                auto indices2(indices<P2, Dim>());
                Index_t<P2> index2(indices2[j]);

                Index_t<P1 + P2> index3;

                size_t m = 0;
                for (; m < P1; ++m) {
                    index3[m] = index1[m];
                }
                for (size_t n = 0; n < P2; ++m, ++n) {
                    index3[m] = index2[n];
                }
                std::sort(index3.begin(), index3.end());
                retval.at(i).at(j) =
                    index_to_canonical<P1 + P2, Dim>().at(index3);
            }
        }
        return retval;
    }

    template <unsigned int P1, unsigned int P2>
    unsigned int f(unsigned int i, unsigned j)
    {
        static std::array<
            std::array<unsigned int, Trigon<double, P2, Dim>::count>,
            Trigon<double, P1, Dim>::count>
            mapping = calculate_f<P1, P2>();
        return mapping[i][j];
    }

    template <unsigned int New_power, typename Mult_trigon_t, typename Array_t>
    void collect_products(Mult_trigon_t const& t, Array_t& new_terms)
    {
        // this x right = new
        constexpr unsigned int right_power = New_power - Power;
        if (right_power <= t.power()) {
            auto& right_terms(t.template get_subpower<right_power>().terms);
            for (size_t i = 0; i < terms.size(); ++i) {
                for (size_t j = 0; j < right_terms.size(); ++j) {
                    size_t k = f<Power, right_power>(i, j);
                    new_terms[k] += terms[i] * right_terms[j];
                }
            }
        }
        lower.template collect_products<New_power>(t, new_terms);
    }

    Trigon<T, Power, Dim> operator*=(Trigon<T, Power, Dim> const& t)
    {
        if (Power > 1) {
            Terms_t new_terms;
            new_terms.fill(0);
            collect_products<Power>(t, new_terms);
            lower *= t.lower;
            terms = new_terms;
        } else {
            const T this_value = value();
            const T right_value = t.value();
            for(size_t i = 0; i < terms.size(); ++i) {
                terms[i] *= right_value;
                terms[i] += t.terms[i] * this_value;
            }
            value() *= right_value;
        }

        return *this;
    }

    Trigon<T, Power, Dim> operator*=(T val)
    {
        for (auto&& c : terms) {
            c *= val;
        }
        lower *= val;
        return *this;
    }

    Trigon<T, Power, Dim> operator*(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval *= t;
        return retval;
    }

    Trigon<T, Power, Dim> operator*(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval *= val;
        return retval;
    }

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
            for (size_t i = 0; i < terms.size(); ++i) {
                terms[i] *= right_value;
                terms[i] -= t.terms[i] * this_value;
                terms[i] *= inv_right_value2;
            }
            value() /= right_value;
        }
        return *this;
    }

    Trigon<T, Power, Dim> operator/=(T val)
    {
        for (auto&& c : terms) {
            c /= val;
        }
        lower /= val;
        return *this;
    }

    Trigon<T, Power, Dim> operator/(Trigon<T, Power, Dim> const& t) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval /= t;
        return retval;
    }

    Trigon<T, Power, Dim> operator/(T val) const
    {
        Trigon<T, Power, Dim> retval(*this);
        retval /= val;
        return retval;
    }
};

template <typename T, unsigned int Dim>
class Trigon<T, 0, Dim>
{
public:
    static constexpr unsigned int count = 1;
    typedef std::array<T, count> Terms_t;
    Terms_t terms;

    Trigon() { terms[0] = 0; }

    Trigon(T val) { terms[0] = val; }

    static constexpr unsigned int power() { return 0; }

    const T value() const { return terms[0]; }

    template <unsigned int Subpower>
    Trigon<T, Subpower, Dim>& get_subpower()
    {
        return *this;
    }

    template <unsigned int Subpower>
    const Trigon<T, Subpower, Dim>& get_subpower() const
    {
        return *this;
    }

    Trigon<T, 0, Dim>& operator+=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] += t.terms[0];
        return *this;
    }

    Trigon<T, 0, Dim>& operator+=(T t)
    {
        terms[0] += t;
        return *this;
    }

    Trigon<T, 0, Dim> operator+(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval += t;
        return retval;
    }

    Trigon<T, 0, Dim> operator+(T t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval += t;
        return retval;
    }

    Trigon<T, 0, Dim> operator-() const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval.terms[0] = -terms[0];
        return retval;
    }

    Trigon<T, 0, Dim>& operator-=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] -= t.terms[0];
        return *this;
    }

    Trigon<T, 0, Dim>& operator-=(T t)
    {
        terms[0] -= t;
        return *this;
    }

    Trigon<T, 0, Dim> operator-(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval -= t;
        return retval;
    }

    Trigon<T, 0, Dim> operator-(T t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval -= t;
        return retval;
    }

    template <unsigned int New_power, typename Mult_trigon_t, typename Array_t>
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

    Trigon<T, 0, Dim> operator*=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] *= t.terms[0];
        return *this;
    }

    Trigon<T, 0, Dim> operator*=(T val)
    {
        terms[0] *= val;
        return *this;
    }

    Trigon<T, 0, Dim> operator*(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval *= t;
        return retval;
    }

    Trigon<T, 0, Dim> operator*(T val) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval *= val;
        return retval;
    }

    Trigon<T, 0, Dim> operator/=(Trigon<T, 0, Dim> const& t)
    {
        terms[0] /= t.terms[0];
        return *this;
    }

    Trigon<T, 0, Dim> operator/=(T val)
    {
        terms[0] /= val;
        return *this;
    }

    Trigon<T, 0, Dim> operator/(Trigon<T, 0, Dim> const& t) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval /= t;
        return retval;
    }

    Trigon<T, 0, Dim> operator/(T val) const
    {
        Trigon<T, 0, Dim> retval(*this);
        retval /= val;
        return retval;
    }
};

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<T, Power, Dim>
operator+(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval(t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
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
Trigon<std::complex<T>, Power, Dim>
operator+(Trigon<std::complex<T>, Power, Dim> const& t1,
          Trigon<T, Power, Dim> const& t2)
{
    return t2 + t1;
}

template <typename T, unsigned int Dim>
Trigon<T, 0, Dim>
operator+(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval(t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Dim>
Trigon<std::complex<T>, 0, Dim> operator+(
    Trigon<T, 0, Dim> const& t1, Trigon<std::complex<T>, 0, Dim> const& t2)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = t1.terms[0] + t2.terms[0];
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<T, Power, Dim>
operator-(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval(-t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Dim>
Trigon<T, 0, Dim>
operator-(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval(-t);
    retval += val;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
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
Trigon<std::complex<T>, 0, Dim> operator-(
    Trigon<T, 0, Dim> const& t1, Trigon<std::complex<T>, 0, Dim> const& t2)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = t1.terms[0] - t2.terms[0];
    return retval;
}

template <typename T, unsigned int Dim>
Trigon<std::complex<T>, 0, Dim> operator-(
    Trigon<std::complex<T>, 0, Dim> const& t1, Trigon<T, 0, Dim> const& t2)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = t1.terms[0] - t2.terms[0];
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<T, Power, Dim> operator*(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval(t);
    retval *= val;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<std::complex<T>, Power, Dim> operator*(std::complex<T> val,
                                              Trigon<T, Power, Dim> const& t)
{
    Trigon<std::complex<T>, Power, Dim> retval;
    for (size_t i = 0; i < retval.terms.size(); ++i) {
        retval.terms[i] = val * t.terms[i];
    }
    retval.lower = val * t.lower;
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<std::complex<T>, Power, Dim> operator*(
    Trigon<std::complex<T>, Power, Dim> const& t1,
    Trigon<T, Power, Dim> const& t2)
{
    return t1 * (std::complex<T>(1.0, 0.0) * t2);
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<std::complex<T>, Power, Dim> operator*(
    Trigon<T, Power, Dim> const& t1,
    Trigon<std::complex<T>, Power, Dim> const& t2)
{
    return (std::complex<T>(1.0, 0.0) * t1) * t2;
}

template <typename T, unsigned int Dim>
Trigon<T, 0, Dim> operator*(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval(t);
    retval *= val;
    return retval;
}

template <typename T, unsigned int Dim>
Trigon<std::complex<T>, 0, Dim> operator*(std::complex<T> val,
                                          Trigon<T, 0, Dim> const& t)
{
    Trigon<std::complex<T>, 0, Dim> retval;
    retval.terms[0] = val * t.terms[0];
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<T, Power, Dim>
operator/(T val, Trigon<T, Power, Dim> const& t)
{
    Trigon<T, Power, Dim> retval;
    retval.template get_subpower<0>().terms[0] = val;
    retval /= t;
    return retval;
}

template <typename T, unsigned int Dim>
Trigon<T, 0, Dim>
operator/(T val, Trigon<T, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval;
    retval.terms[0] = val / t.terms[0];
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
typename std::enable_if<(Power == 1), void>::type
calculate_partial(Trigon<T, Power, Dim> const& source, size_t index,
                  Trigon<T, Power - 1, Dim>& dest)
{
    dest.terms[0] = source.terms[index];
}

template <typename T, unsigned int Power, unsigned int Dim>
typename std::enable_if<(Power > 1), void>::type
calculate_partial(Trigon<T, Power, Dim> const& source, size_t index,
                  Trigon<T, Power - 1, Dim>& dest)
{
    auto source_indices(indices<Power, Dim>());
    for (size_t i = 0; i < source_indices.size(); ++i) {
        auto& source_index(source_indices[i]);
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
            dest.terms[index_to_canonical<Power - 1, Dim>().at(dest_index)] =
                exponent * source.terms[i];
        }
    }
    calculate_partial<T, Power - 1, Dim>(source.lower, index, dest.lower);
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<T, Power - 1, Dim>
partial_deriv(Trigon<T, Power, Dim> const& t, size_t index)
{
    Trigon<T, Power - 1, Dim> retval;
    calculate_partial<T, Power, Dim>(t, index, retval);
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
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
Trigon<T, 0, Dim> real(Trigon<std::complex<T>, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval;
    retval.terms[0] = t.terms[0].real();
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
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
Trigon<T, 0, Dim> imag(Trigon<std::complex<T>, 0, Dim> const& t)
{
    Trigon<T, 0, Dim> retval;
    retval.terms[0] = t.terms[0].imag();
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
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
using Derivatives_t = std::array<T, max_power + 1>;

template <typename T, unsigned int Power, unsigned int Dim>
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
Trigon<T, Power, Dim>
sin(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, sin_derivatives(t.value(), Power));
}

template <typename T>
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
Trigon<T, Power, Dim>
cos(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, cos_derivatives(t.value(), Power));
}

template <typename T>
T
qpow(T x, int i)
{
    T retval = 1;
    while (i--) {
        retval *= x;
    }
    return retval;
}

template <typename T>
Derivatives_t<T>
tan_derivatives(T x, unsigned int power)
{
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
Derivatives_t<T>
sqrt_derivatives(T x, unsigned int power)
{
    Derivatives_t<T> retval;
    T invsqrtx;
    retval[0] = std::sqrt(x);
    if (power > 0) {
        invsqrtx = 1 / retval[0];
        retval[1] = 0.5 * invsqrtx;
    }
    if (power > 1) {
        retval[2] = -0.5 * invsqrtx * retval[1];
    }
    if (power > 2) {
        retval[3] = -1.5 * invsqrtx * retval[2];
    }
    if (power > 3) {
        retval[4] = -2.5 * invsqrtx * retval[3];
    }
    if (power > 4) {
        retval[5] = -3.5 * invsqrtx * retval[4];
    }
    if (power > 5) {
        retval[6] = -4.5 * invsqrtx * retval[5];
    }
    if (power > 6) {
        retval[7] = -5.5 * invsqrtx * retval[6];
    }
    return retval;
}

template <typename T, unsigned int Power, unsigned int Dim>
Trigon<T, Power, Dim>
sqrt(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, sqrt_derivatives(t.value(), Power));
}

template <typename T>
Derivatives_t<T>
asin_derivatives(T x, unsigned int power)
{
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
Trigon<T, Power, Dim>
asin(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, asin_derivatives(t.value(), Power));
}

template <typename T>
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
Trigon<T, Power, Dim>
log(Trigon<T, Power, Dim> const& t)
{
    return generic_transcendental(t, log_derivatives(t.value(), Power));
}

#endif // TRIGON_H
