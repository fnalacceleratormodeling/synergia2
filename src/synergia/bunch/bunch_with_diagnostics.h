#ifndef BUNCH_DIAG_H_
#define BUNCH_DIAG_H_

#include <mpi.h>
#include "boost/shared_ptr.hpp"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/multi_diagnostics.h"

class Bunch_with_diagnostics
{ 
private:
    Bunch_sptr bunch_sptr; 
    Multi_diagnostics  per_step_diagnostics;
    Multi_diagnostics  per_turn_diagnostics;
    
    void
    check_bunch_pointer();
    
    void
    construct(Multi_diagnostics diagnostics_step,
                         Multi_diagnostics diagnostics_turn);
public:    
     Bunch_with_diagnostics(Bunch_sptr bunch_sptr, Diagnostics_sptr diagnostics_step_sptr, 
                        Diagnostics_sptr diagnostics_turn_sptr); 
                        
     Bunch_with_diagnostics(Bunch_sptr bunch_sptr, Multi_diagnostics step_diagnostics,
                        Multi_diagnostics  turn_diagnostics);
                        
    
         
     Bunch_sptr 
     get_bunch_sptr();
      
     Multi_diagnostics &
     get_per_step_diagnostics();
     
     Multi_diagnostics &
     get_per_turn_diagnostics();
     
     ~Bunch_with_diagnostics();  
};

typedef boost::shared_ptr<Bunch_with_diagnostics > Bunch_with_diagnostics_sptr;

#endif /* BUNCH_DIAG_H_ */
