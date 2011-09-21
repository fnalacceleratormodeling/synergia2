#include "bunch_with_diagnostics.h"


void
Bunch_with_diagnostics::construct(Multi_diagnostics multi_diagnostics_step,
                         Multi_diagnostics multi_diagnostics_turn)
{   

    this->per_step_diagnostics=multi_diagnostics_step;
    this->per_turn_diagnostics=multi_diagnostics_turn;                 
    check_bunch_pointer();
}






Bunch_with_diagnostics::Bunch_with_diagnostics(Bunch_sptr bunch_sptr, Multi_diagnostics multi_diagnostics_step,
                         Multi_diagnostics multi_diagnostics_turn): 
                         bunch_sptr(bunch_sptr)
                         
{   
        construct(multi_diagnostics_step,multi_diagnostics_turn);
}


Bunch_with_diagnostics::Bunch_with_diagnostics(Bunch_sptr bunch_sptr, 
                            Diagnostics_sptr diagnostics_step_sptr, 
                            Diagnostics_sptr diagnostics_turn_sptr):
                            bunch_sptr(bunch_sptr)

{    
        Multi_diagnostics multi_diagnostics_step;
        multi_diagnostics_step.append(diagnostics_step_sptr);
        Multi_diagnostics multi_diagnostics_turn;
        multi_diagnostics_turn.append(diagnostics_turn_sptr); 
        construct(multi_diagnostics_step,multi_diagnostics_turn);
 }
 


 
void 
Bunch_with_diagnostics::check_bunch_pointer()
{
        for (Multi_diagnostics::iterator dit = this->per_step_diagnostics.begin(); dit
                !=  this->per_step_diagnostics.end(); ++dit) {
            if ((*dit)->get_bunch_sptr() != bunch_sptr) 
                throw std::runtime_error(
                 "the per_step_diagnostics bunch pointer and the bunch are different objects");

        }

       for (Multi_diagnostics::iterator dit = this->per_turn_diagnostics.begin(); dit
                !=  this->per_turn_diagnostics.end(); ++dit) {
            if ((*dit)->get_bunch_sptr() != bunch_sptr) 
                throw std::runtime_error(
                 "the per_turn_diagnostics bunch pointer and the bunch are different objects");

        }

} 
 
 
Bunch_sptr 
Bunch_with_diagnostics::get_bunch_sptr()
{ 
return bunch_sptr;
}

   
Multi_diagnostics &
Bunch_with_diagnostics::get_per_step_diagnostics()
{
    return per_step_diagnostics;
}
   
   
Multi_diagnostics &
Bunch_with_diagnostics::get_per_turn_diagnostics()
{
    return per_turn_diagnostics;
}
             
Bunch_with_diagnostics:: ~Bunch_with_diagnostics()
 {
 }
      
      
/* */
