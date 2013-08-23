#include <iostream>
#include "step.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/collective/impedance.h"
#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/get_pointer.hpp>
#include "synergia/collective/impedance.h"
#include "synergia/bunch/period.h"

Step::Step(double length) :
    operators(), time_fractions(), length(length)
{
}

Step::Step()
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
Step::apply(Bunch & bunch, int verbosity,
        Diagnosticss const& per_operator_diagnostics,
        Diagnosticss const& per_operation_diagnostics, Logger & logger)
{
    double t_total = simple_timer_current();
    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (Operators::const_iterator it = operators.begin();
            it != operators.end(); ++it) {
        // time [s] in accelerator frame
        double time = length / (bunch.get_reference_particle().get_beta()
                * pconstants::c);        
        double t0 = MPI_Wtime();
        double t = simple_timer_current();
        (*it)->apply(bunch, (*fractions_it) * time, *this, verbosity,
                per_operation_diagnostics, logger);
        std::string label("step_apply-" + (*it)->get_type() + "_operator_apply");
        t = simple_timer_show(t, label.c_str());
        double t1 = MPI_Wtime();
        if (verbosity > 2) {
            logger << "Step: operator: name = " << (*it)->get_name()
                    << ", type = " << (*it)->get_type() << ", time = "
                    << std::fixed << std::setprecision(3) << t1 - t0 << "s_n"
                    << std::endl;
        }

        t = simple_timer_current();
        for (Diagnosticss::const_iterator itd =
                per_operator_diagnostics.begin();
                itd != per_operator_diagnostics.end(); ++itd) {
            (*itd)->update_and_write();
        }
        t = simple_timer_show(t, "diagnostics-operator");

        if (bunch.is_z_periodic()) {
            double plength = bunch.get_z_period_length();
            apply_longitudinal_periodicity(bunch, plength);
        }
        ++fractions_it;
    }
    t_total = simple_timer_show(t_total, "step_apply-total");
}

void
Step::apply(Bunch_train & bunch_train, int verbosity,
        Train_diagnosticss const& per_operator_train_diagnosticss,
        Train_diagnosticss const& per_operation_train_diagnosticss, Logger & logger)
{
    // time [s] in accelerator frame
    double time = length
            / (bunch_train.get_bunches()[0]->get_reference_particle().get_beta()
                    * pconstants::c);
    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (Operators::const_iterator it = operators.begin();
            it != operators.end(); ++it) {
        double t0 = MPI_Wtime();
        (*it)->apply(bunch_train, (*fractions_it) * time, *this, verbosity,
                per_operation_train_diagnosticss, logger);
        double t1 = MPI_Wtime();
        if (verbosity > 2) {
            logger << "Step: operator: name = " << (*it)->get_name()
                    << ", type = " << (*it)->get_type() << ", time = "
                    << std::fixed << std::setprecision(3) << t1 - t0 << "s_n"
                    << std::endl;
        }

        double t = simple_timer_current();
        size_t num_bunches = bunch_train.get_size();
        for (int i = 0; i < num_bunches; ++i) {
            for (Diagnosticss::const_iterator itd =
                    per_operator_train_diagnosticss.at(i).begin();
                    itd != per_operator_train_diagnosticss.at(i).end(); ++itd) {
                (*itd)->update_and_write();
            }
        }
        t = simple_timer_show(t, "diagnostics-operator");
        // jfa: what should we do here? Move particles between bunches?
	// am: this should be changed soon, for now assume the effect is small
        for (int i = 0; i < num_bunches; ++i) {
	   if (bunch_train.get_bunches().at(i)->get_comm().has_this_rank()) {
                Bunch_sptr bunch_sptr=bunch_train.get_bunches().at(i);
		double plength=bunch_sptr->get_z_period_length();
                if (bunch_sptr->is_z_periodic()){              
                    apply_longitudinal_periodicity(*bunch_sptr, plength);
                } 
                else{                 
                    apply_zcut(*bunch_sptr, plength);
                }
	   }                  
        }           
        ++fractions_it;
    }
}



Operators const&
Step::get_operators() const
{
    return operators;
}

Operators &
Step::get_operators()
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



void
Step::print(int index) const
{
    std::cout << "step " << index << ":\n";
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        (*it)->print();
    }
}

template<class Archive>
    void
    Step::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(operators);
        ar & BOOST_SERIALIZATION_NVP(time_fractions);
        ar & BOOST_SERIALIZATION_NVP(length);       
    }

template
void
Step::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Step::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Step::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Step::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
