#ifndef LATTICE_ELEMENT_SLICE_H_
#define LATTICE_ELEMENT_SLICE_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "components/lattice/lattice_element.h"

class Lattice_element_slice
{
private:
    Lattice_element * element_ptr;
    bool whole;
    bool left_edge;
    bool right_edge;
    double left;
    double right;

public:
    Lattice_element_slice(Lattice_element & lattice_element);
    Lattice_element_slice(Lattice_element & lattice_element, double left,
            double right);
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
    void
    print() const;
};

typedef boost::shared_ptr<Lattice_element_slice > Lattice_element_slice_sptr;
typedef std::list<Lattice_element_slice_sptr > Lattice_element_slices;

#endif /* LATTICE_ELEMENT_SLICE_H_ */
