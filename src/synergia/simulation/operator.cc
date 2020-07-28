#include <iostream>
#include "synergia/utils/simple_timer.h"
#include "operator.h"
//#include "aperture_operation.h"

Operator::Operator(std::string const& name, std::string const& type) :
        name(name), type(type)
{
}

Operator::Operator()
{
}

std::string const&
Operator::get_name() const
{
    return name;
}

std::string const&
Operator::get_type() const
{
    return type;
}

void
Operator::apply(Bunch_train & bunch_train, double time_step, Step & step,
        int verbosity, Train_diagnosticss const& per_operation_diagnosticss,
        Logger & logger)
{
    Bunches bunches(bunch_train.get_bunches());
    size_t num_bunches = bunch_train.get_size();
    for (std::size_t i = 0; i < num_bunches; ++i)
        if (bunches.at(i)->get_comm().has_this_rank()) {
            apply(*bunches.at(i), time_step, step, verbosity,
                    per_operation_diagnosticss.at(i), logger);
        }
}


void
Operator::apply(Bunch_train & bunch_train, double time_step, Step & step, int verbosity,
          Train_diagnosticss const& per_operation_diagnosticss, 
          Propagate_actions * propagate_actions_ptr, Stepper & stepper, int step_count,  int turn, 
          Logger & logger)
{
    Bunches bunches(bunch_train.get_bunches());
    size_t num_bunches = bunch_train.get_size();
    for (std::size_t i = 0; i < num_bunches; ++i)
        if (bunches.at(i)->get_comm().has_this_rank()) {  
            propagate_actions_ptr->operator_begin_action(stepper, step, *this, step_count,  turn,  i);             
            apply(*bunches.at(i), time_step, step, verbosity,
                    per_operation_diagnosticss.at(i), logger);                 
        }
}

void
Operator::print() const
{
    std::cout << type << " operator: " << name << std::endl;
}



#if 0
void
Operator::apply_train(Bunch_with_diagnostics_train & bunch_diag_train,
        double time_step, Step & step)
{
    for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
        if (bunch_diag_train.is_on_this_rank(index)) {
            Bunch_sptr
            bunch_sptr =
            bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr();
            apply(*bunch_sptr, time_step, step);
        }
    }
}
#endif

template<class Archive>
    void
    Operator::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(name);
        ar & BOOST_SERIALIZATION_NVP(type);
    }

template
void
Operator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Operator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Operator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Operator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Operator::~Operator()
{
}


