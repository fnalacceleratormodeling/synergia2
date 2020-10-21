#ifndef DUMMY_COLLECTIVE_OPERATOR_H
#define DUMMY_COLLECTIVE_OPERATOR_H

#include "synergia/simulation/operator.h"
#include "synergia/simulation/collective_operator_options.h"

class Dummy_collective_operator;

// Dummy collective
struct Dummy_CO_options
    : public CO_base_options<
        Dummy_CO_options, 
        Dummy_collective_operator>
{
    template<class Archive>
    void serialize(Archive & ar)
    { ar(cereal::base_class<CO_base_options>(this)); }
};

CEREAL_REGISTER_TYPE(Dummy_CO_options)

class Dummy_collective_operator : public Collective_operator
{
private:

    virtual void apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger) override
    { }

public:

    Dummy_collective_operator(Dummy_CO_options const & ops)
        : Collective_operator("dummy collective", 1.0)
    { }
};

#endif
