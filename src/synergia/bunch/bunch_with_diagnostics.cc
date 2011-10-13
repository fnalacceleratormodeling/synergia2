#include "bunch_with_diagnostics.h"

Bunch_with_diagnostics::Bunch_with_diagnostics(Bunch_sptr bunch_sptr,  Propagate_actions_sptr diagnostics_actions_sptr):
 bunch_sptr(bunch_sptr), diagnostics_actions_sptr(diagnostics_actions_sptr)
{  
}


 
Bunch_sptr 
Bunch_with_diagnostics::get_bunch_sptr()
{ 
return bunch_sptr;
}


Propagate_actions_sptr &
Bunch_with_diagnostics::get_diagnostics_actions_sptr()
{
    return diagnostics_actions_sptr;
}           
             

Bunch_with_diagnostics:: ~Bunch_with_diagnostics()
 {
 }
      
      
/* */
