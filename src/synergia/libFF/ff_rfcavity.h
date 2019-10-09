#ifndef FF_RFCAVITY_H
#define FF_RFCAVITY_H

#include "ff_element.h"

class FF_rfcavity : public FF_element
{
public:
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_RFCAVITY_H
