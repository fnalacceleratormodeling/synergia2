#ifndef FF_MULTIPOLE_H
#define FF_MULTIPOLE_H

#include "ff_element.h"

class FF_multipole : public FF_element
{
public:
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_MULTIPOLE_H
