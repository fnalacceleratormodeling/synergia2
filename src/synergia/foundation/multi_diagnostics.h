#ifndef MULTI_DIAGNOSTICS_H_
#define MULTI_DIAGNOSTICS_H_
#include "synergia/foundation/generalized_diagnostics.h"
#include <string>
#include <list>

/// Multi_diagnostics contains a list of Generalized_diagnostics_sptrs
class Multi_diagnostics
{
private:
    std::list<Generalized_diagnostics_sptr > diagnostics;
public:
    /// Construct an empty list of Generalized_diagnostics_sptrs
    Multi_diagnostics();

    /// Append a Generalized_diagnostics_sptr to the list
    /// @param diagnostics_sptr the Generalized_diagnostics_sptr
    void
    append(Generalized_diagnostics_sptr diagnostics_sptr);

    /// The same as append -- included for notational consistency with C++
    /// @param diagnostics_sptr the Generalized_diagnostics_sptr
    void
    push_back(Generalized_diagnostics_sptr diagnostics_sptr);

    /// A convenience definition of the list iterator. Not relevant for
    /// Python.
    typedef std::list<Generalized_diagnostics_sptr >::iterator iterator;

    /// A convenience definition of the list iterator begin. Not relevant for
    /// Python.
    iterator
    begin();

    /// A convenience definition of the list iterator end. Not relevant for
    /// Python.
    iterator
    end();
};

///// A pre-defined dummy Multi_diagnostics
Multi_diagnostics
no_diagnostics();

#endif /* MULTI_DIAGNOSTICS_H_ */
