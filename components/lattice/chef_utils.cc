#include "chef_utils.h"
#include <iostream>

void
print_chef_beamline(BmlPtr beamline_sptr)
{
    for (beamline::const_iterator it = beamline_sptr->begin(); it
            != beamline_sptr->end(); ++it) {
        std::cout << (*it)->Name() << "(" << (*it)->Type() << "): Length="
                << (*it)->Length() << ", Strength=" << (*it)->Strength()
                << std::endl;
    }
}
