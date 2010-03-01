#ifndef LATTICE_ELEMENT_SLICE_H_
#define LATTICE_ELEMENT_SLICE_H_

#include "components/lattice/lattice_element.h"
#include <list>

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
};

#endif /* LATTICE_ELEMENT_SLICE_H_ */
