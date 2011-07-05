#include <iostream>
#include "step.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/collective/impedance.h"
#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/get_pointer.hpp>
#include "synergia/collective/impedance.h"


Step::Step(double length) :
    operators(), time_fractions(), length(length)
{

}

void
Step::append(Operator_sptr operator_sptr, double time_fraction)
{
    operators.push_back(operator_sptr);
    time_fractions.push_back(time_fraction);
}

void
Step::append(Operators const& the_operators, double time_fraction)
{
    Operators tmp(the_operators);
    operators.splice(operators.end(), tmp);
    for (Operators::const_iterator it = the_operators.begin(); it
            != the_operators.end(); ++it) {
        time_fractions.push_back(time_fraction);
    }
}

void
Step::apply(Bunch & bunch)
{ 
    
    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        // time [s] in accelerator frame
        double time = length / (bunch.get_reference_particle().get_beta()
                * pconstants::c);
        if ((*it)->get_name()=="impedance") {
            MArray1d bunch_means=Diagnostics::calculate_mean(bunch);
            Bunch_means bi;
            bi.x_mean=bunch_means[0];
            bi.y_mean=bunch_means[2];
            bi.z_mean=bunch_means[4];
            bi.n_part=bunch.get_total_num();
            stored_bunches.push_front(bi);

            int nstored=(reinterpret_cast<Impedance*>(boost::get_pointer(*it)))->get_nstored_turns(); 
            if (stored_bunches.size()>nstored) stored_bunches.pop_back();
          //  std::cout<<"name ="<< (*it)->get_name()<<" stored dim "<<stored_bunches.size()<<std::endl; 
           
         }
        (*it)->apply(bunch, (*fractions_it) * time, *this);
        ++fractions_it;
    }
}

Operators const&
Step::get_operators() const
{
    return operators;
}

std::list<double > const&
Step::get_time_fractions() const
{
    return time_fractions;
}

double
Step::get_length() const
{
    return length;
}

std::list<Bunch_means>  Step::get_stored_bunches() const
{
        return stored_bunches;
 }

void
Step::print(int index) const
{
    std::cout << "step " << index << ":\n";
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        (*it)->print();
    }
}
