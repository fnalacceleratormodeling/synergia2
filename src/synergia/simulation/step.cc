#include <iostream>

#include "step.h"
//#include "synergia/utils/simple_timer.h"
#include "synergia/foundation/physical_constants.h"

Step::Step(double length) 
: operators()
, step_betas()
, length(length)
{
}

#if 0
Step::Step()
: operators()
, step_betas()
, length(0.0)
{
}
#endif

void Step::apply(Bunch_simulator & simulator, Logger & logger) const
{
    // time [s] in accelerator frame
    double ref_beta = simulator.get_bunch(0, 0).get_reference_particle().get_beta();
    double time = length / (ref_beta * pconstants::c);

    for (auto const & op : operators)
    {
        auto const & opr = op.first;
        double fraction = op.second;

        double t0 = MPI_Wtime();

        // operator apply
        opr->apply(simulator, fraction * time, logger);

        double t1 = MPI_Wtime();

        //if (verbosity > 2) 
        {
            logger << "Step: operator: name = " << opr->get_name()
                    << ", type = " << opr->get_type() << ", time = "
                    << std::fixed << std::setprecision(3) << t1 - t0 << "s"
                    << std::endl;
        }

        //double t = simple_timer_current();

        // per operator diagnostics action
        simulator.diag_action_operator(*opr);

        //t = simple_timer_show(t, "diagnostics-operator");      
    }
}

double
Step::get_length() const
{
    return length;
}

#if 0
void
Step::set_betas(double betax, double betay)
{
 this->step_betas.push_back(betax);
 this->step_betas.push_back(betay);
}

std::vector<double >
Step::get_betas()
{
 return step_betas;
}

void
Step::print(int index) const
{
#if 0
    std::cout << "step " << index << ":\n";
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        (*it)->print();
    }
#endif
}

#if 0
template<class Archive>
    void
    Step::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(operators);
        ar & BOOST_SERIALIZATION_NVP(time_fractions);
        ar & BOOST_SERIALIZATION_NVP(length); 
        ar & BOOST_SERIALIZATION_NVP(step_betas);
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
#endif
#endif
