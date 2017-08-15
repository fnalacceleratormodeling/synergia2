#ifndef FAST_MAPPING_H_
#define FAST_MAPPING_H_

#include "synergia/foundation/reference_particle.h"
#include "synergia/bunch/bunch.h"
#include "synergia/lattice/chef_utils.h"

#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#include "mxyzptlk/Mapping.h"
#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif

#include <vector>
#include <list>
#include <fstream>
#include <string>

#include "synergia/utils/serialization.h"

// The interface in this comment is out of date
// A Fast_mapping_term represents one term in a polynomial expansion of order
// "order". It contains a coefficient "coeff" and a vector "i" of length
// "order" containing the vector indices of the dependent variable.
//
// Examples:
//   3.14 * p_1
//   => Fast_mapping_term fmt(1); fmt.coeff() = 3.14; fmt.index(0) = 1;
//
//   1.7724 * p_0^3 (i.e., p_0 cubed)
//   => Fast_mapping_term fmt(3); fmt.coeff() = 1.7724;
//          fmt.index(0) = 0; fmt.index(1) = 0; fmt.index(2) = 0;
//
//   2.2 * p_0*p_2*p_4
//   => Fast_mapping_term fmt(3); fmt.coeff() = 2.2;
//          fmt.index(0) = 0; fmt.index(1) = 2; fmt.index(2) = 4;
//
template < typename T> 
class Fast_mapping_term
{
private:
    T the_coeff;
    std::vector<int > i;
    int the_order;
public:
    Fast_mapping_term(int order);
    Fast_mapping_term(std::ifstream & stream);
    /// Default constructor for serialization use only
    Fast_mapping_term();
    Fast_mapping_term(Fast_mapping_term const& fast_mapping_term);
    inline int
    order() const
    {
        return the_order;
    }
    ;
    inline T &
    coeff()
    {
        return the_coeff;
    }
    ;
    inline T const&
    coeff() const
    {
        return the_coeff;
    }
    ;
    inline int &
    index(int which)
    {
        return i[which];
    }
    ;
    inline int const&
    index(int which) const
    {
        return i[which];
    }
    ;
    void
    write_to_stream(std::ostream & stream) const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~Fast_mapping_term();
};

//typedef std::list<Fast_mapping_term > Fast_mapping_terms; // syndoc:include


// The interface in this comment is out of date
// Fast_mapping is a sparse container for a collection of arbitrary-order
// terms in a polynomial expansion taking the six-dimensional phase space
// variable g to a new variable f:
//
// f^{i} = C_{0}^{i} + C_{1}^{ij} g^{j} +
//         C_{2}^{ijk} g^{j}g^{k} +
//         C_{3}^{ijkl}g^{j}g^{k}g^{l} + ...
//
// Only the non-zero terms are stored.
// The "terms" member contains a vector of vectors of lists of
// "Fast_mapping_terms", each of which corresponds to one term in the
// fully expanded version of the above polynomial. The details of the
// storage are designed for optimized application of the mappings, not
// not necessarily for clarity. The structure is described below.
//
// The outermost vector in terms runs from 0 to 5, corresponding to the
// components x px y py t pt of the output vector, $i$ in the equation
// above.
// For each $i$, There is a vector running from 0 to N+1, where N is the
// highest order term present. Each component of this vector contains
// a list of all the non-zero terms at that order. The Fast_mapping_term's
// themselves are defined near the definition of the class.
//
// Example, using only a single value of $i$:
//
// f^1 = 1.1 + 2.2*g^1 - 3.0*g^2 + 4.0*g^4*g^5
//
// Fast_mapping_term tmp_term0(0);
// tmp_term0.coeff = 1.1;
// terms.at(1).at(0).push_back(tmp_term0);
//
// Fast_mapping_term tmp_term1(1);
// tmp_term1.coeff = 2.2;
// tmp_term1.i[0] = 1;
// terms.at(1).at(1).push_back(tmp_term1);
// tmp_term1.coeff = 3.0;
// tmp_term1.i[0] = 2;
// terms.at(1).at(1).push_back(tmp_term1);
//
// Fast_mapping_term tmp_term2(2);
// tmp_term2.coeff = 4.0;
// tmp_term2.i[0] = 4;
// tmp_term2.i[1] = 5;
// terms.at(1).at(2).push_back(tmp_term2);
//

template < typename T> 
class TFast_mapping
{
private:
    std::vector<std::vector<std::list<Fast_mapping_term<T>  > > > terms;
    int order;
    double length;
    void
    init(int order);
public:
    TFast_mapping(int order);
    TFast_mapping(std::string const& filename);
    TFast_mapping(TMapping<T> const& chef_mapping, double mapping_length);
    /// Default constructor for serialization use only
    TFast_mapping();
    void
    set_length(double length);
    double
    get_length() const;
    int
    get_order() const;
    void
    add_term(int index, Fast_mapping_term<T> const& term);
    std::vector<std::vector<std::list<Fast_mapping_term<T>  > > > const&
    get_terms() const;
    void
    apply(Bunch & bunch);
    void
    apply(boost::multi_array_ref<T, 1 > coords);   
//     T
//     apply(int icoord,  boost::multi_array_ref<T, 1 > coords);
    
    std::string
    as_string() const;
    void
    write_to_file(std::string const& filename);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

typedef TFast_mapping<double>  Fast_mapping;
typedef  boost::shared_ptr<Fast_mapping > Fast_mapping_sptr;
typedef TFast_mapping<std::complex<double > >  CFast_mapping;
typedef  boost::shared_ptr<CFast_mapping > CFast_mapping_sptr;

#include "synergia/simulation/fast_mapping.tcc"
#endif /* FAST_MAPPING_H_ */
