#ifndef DUMMY_COLLECTIVE_OPERATOR_H
#define DUMMY_COLLECTIVE_OPERATOR_H

#include "synergia/simulation/operator.h"
#include "synergia/simulation/collective_operator_options.h"

// Dummy collective
struct Dummy_CO_options : public CO_options
{
    Collective_operator * create_operator() const override;
    CO_options * clone() const override { return new Dummy_CO_options(*this); }

    template<class Archive>
    void serialize(Archive & ar)
    { ar(cereal::base_class<CO_options>(this)); }
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

inline Collective_operator *
Dummy_CO_options::create_operator() const
{ return new Dummy_collective_operator(*this); }


#endif
