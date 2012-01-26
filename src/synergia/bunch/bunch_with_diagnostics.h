#ifndef BUNCH_DIAG_H_
#define BUNCH_DIAG_H_

#include <mpi.h>
#include "boost/shared_ptr.hpp"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/bunch/diagnostics.h"


class Propagate_actions;
typedef boost::shared_ptr<Propagate_actions > Propagate_actions_sptr;
class Standard_diagnostics_actions;
typedef boost::shared_ptr<Standard_diagnostics_actions > Standard_diagnostics_actions_sptr;

class Bunch_with_diagnostics
{
private:
    Bunch_sptr bunch_sptr;
    Standard_diagnostics_actions_sptr diagnostics_actions_sptr;

public:

     void
     check_bunch_pointer_in_diagnostics() const;

     Bunch_with_diagnostics(Bunch_sptr bunch_sptr,  Standard_diagnostics_actions_sptr diagnostics_actions_sptr);

     void
     add_per_step_diagnostics(Diagnostics_sptr diagnostics_sptr);

     void
     add_per_turn_diagnostics(Diagnostics_sptr diagnostics_sptr);

     Bunch_sptr const
     get_bunch_sptr();

     Standard_diagnostics_actions_sptr
     get_diagnostics_actions_sptr() const;


    /// Get the bunch communicator
     Commxx const&
     get_comm() const;

     ~Bunch_with_diagnostics();
};

typedef boost::shared_ptr<Bunch_with_diagnostics > Bunch_with_diagnostics_sptr;

#endif /* BUNCH_DIAG_H_ */
