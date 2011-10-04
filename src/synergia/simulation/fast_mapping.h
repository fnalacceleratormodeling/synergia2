#ifndef FAST_MAPPING_H_
#define FAST_MAPPING_H_

#include "synergia/foundation/reference_particle.h"
#include "synergia/bunch/bunch.h"
#include "mxyzptlk/Mapping.h"
#include "synergia/lattice/chef_utils.h"

#include <vector>
#include <list>
#include <fstream>
#include <string>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>

// The interface in this comment is out of date
// A Fast_mapping_term represents one term in a polynomial expansion of order
// "order". It contains a coefficient "coeff" and a c-style vector "i" of length
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
class Fast_mapping_term
{
private:
    double the_coeff;
    int *i;
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
    inline double &
    coeff()
    {
        return the_coeff;
    }
    ;
    inline double const&
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
    write_to_stream(std::ofstream & stream) const;
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const
        {
            ar & BOOST_SERIALIZATION_NVP(the_coeff);
            ar & BOOST_SERIALIZATION_NVP(the_order);
            ar & boost::serialization::make_nvp("i",
                    boost::serialization::make_array(i, the_order + 1));
        }
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(the_coeff);
            ar & BOOST_SERIALIZATION_NVP(the_order);
            i = new int[the_order + 1];
            ar & boost::serialization::make_nvp("i",
                    boost::serialization::make_array(i, the_order + 1));
        }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    ~Fast_mapping_term();
};

typedef std::list<Fast_mapping_term > Fast_mapping_terms;

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
class Fast_mapping
{
private:
    std::vector<std::vector<Fast_mapping_terms > > terms;
    int order;
    double length;
    void
    init(int order);
public:
    Fast_mapping(int order);
    Fast_mapping(std::string const& filename);
    Fast_mapping(Reference_particle const& reference_particle,
            Mapping const& chef_mapping, double mapping_length);
    /// Default constructor for serialization use only
    Fast_mapping();
    void
    set_length(double length);
    double
    get_length() const;
    void
    add_term(int index, Fast_mapping_term const& term);
    void
    apply(Bunch & bunch);
    void
    write_to_file(std::string const& filename);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(terms);
            ar & BOOST_SERIALIZATION_NVP(order);
            ar & BOOST_SERIALIZATION_NVP(length);
        }
};

#endif /* FAST_MAPPING_H_ */
