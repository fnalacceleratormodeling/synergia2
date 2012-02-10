#include "multi_diagnostics.h"
#include <sstream>
#include <iomanip>

Multi_diagnostics::Multi_diagnostics() :
    diagnostics()
{
}

void
Multi_diagnostics::append(
        Diagnostics_sptr diagnostics_sptr)
{
    diagnostics.push_back(diagnostics_sptr);
}

void
Multi_diagnostics::push_back(
        Diagnostics_sptr diagnostics_sptr)
{
    diagnostics.push_back(diagnostics_sptr);
}

Multi_diagnostics::iterator
Multi_diagnostics::begin()
{
    return diagnostics.begin();
}

Multi_diagnostics::iterator
Multi_diagnostics::end()
{
    return diagnostics.end();
}

Multi_diagnostics
no_diagnostics()
{
    return Multi_diagnostics();
}
