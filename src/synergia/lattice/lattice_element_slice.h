#ifndef LATTICE_ELEMENT_SLICE_H_
#define LATTICE_ELEMENT_SLICE_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/lattice/lattice_element.h"
#include "synergia/utils/serialization.h"

class Lattice_element_slice
{
private:
    Lattice_element_sptr element_sptr;
    bool whole;
    bool left_edge;
    bool right_edge;
    double left;
    double right;

public:
    Lattice_element_slice(Lattice_element_sptr lattice_element_sptr);
    Lattice_element_slice(Lattice_element_sptr lattice_element_sptr, double left,
            double right);
    // Default constructor for serialization use only
    Lattice_element_slice();
    bool
    is_whole() const;
    bool
    has_left_edge() const;
    bool
    has_right_edge() const;
    double
    get_left() const;
    double
    get_right() const;
    Lattice_element const&
    get_lattice_element() const;
    Lattice_element &
    get_lattice_element();
    std::string
    as_string() const;
    void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

typedef boost::shared_ptr<Lattice_element_slice > Lattice_element_slice_sptr;
typedef std::list<Lattice_element_slice_sptr > Lattice_element_slices;

#endif /* LATTICE_ELEMENT_SLICE_H_ */
